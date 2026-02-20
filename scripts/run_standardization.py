#!/usr/bin/env python
"""
Batch standardization script using the ChEMBL pipeline.

Accepts a single CSV/SDF file or a folder containing CSV and SDF files.
Runs the ChemAudit standardization pipeline on each molecule and writes
results (and failures) to output CSVs.

Usage:
    # Single file — output is a file path
    python scripts/run_standardization.py \
        --input examples/260219_Antimetabolite\ SMILES_H37Rv_Master.csv \
        --output results/standardized_output.csv

    # Folder — output is a directory; each input file gets its own output
    #   e.g. data/mols_A.csv  -> results/mols_A_standardized.csv
    #        data/mols_A.csv  -> results/mols_A_failed.csv  (if any failures)
    #        data/mols_B.sdf  -> results/mols_B_standardized.csv
    python scripts/run_standardization.py \
        --input data/molecules/ \
        --output results/
"""

import argparse
import csv
import glob
import json
import os
import sys
import time

# Add backend to path so we can import the pipeline directly
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "backend"))

from rdkit import Chem

from app.services.standardization.chembl_pipeline import (
    StandardizationOptions,
    standardize_molecule,
)

ID_COL = "id"

FIELDNAMES = [
    ID_COL,
    "source_file",
    "original_smiles",
    "standardized_smiles",
    "success",
    "error",
    "checker_issues",
    "excluded_fragments",
    "stereo_warning",
    "mass_change_pct",
    "formula_original",
    "formula_standardized",
    "steps_summary",
]

FAILED_FIELDNAMES = [
    ID_COL,
    "source_file",
    "original_smiles",
    "error",
    "checker_issues",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ChEMBL standardization pipeline on CSV/SDF files or a folder."
    )
    parser.add_argument(
        "--input", required=True,
        help="Input CSV/SDF file, or a folder containing .csv and .sdf files",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output CSV path (single file) or output directory (folder input with multiple files)",
    )
    parser.add_argument(
        "--smiles-col", default="SMILES",
        help="Default SMILES column name for CSV files (default: SMILES)",
    )
    parser.add_argument(
        "--id-col", default="RegistrationNo",
        help="Default ID column name for CSV files (default: RegistrationNo)",
    )
    parser.add_argument(
        "--column-map", default=None,
        help='Per-file column mapping as JSON: {"file.csv": {"smiles": "COL", "id": "COL"}, ...}',
    )
    parser.add_argument(
        "--include-tautomer", action="store_true",
        help="Enable tautomer canonicalization (WARNING: may lose E/Z stereo)",
    )
    return parser.parse_args()


def get_columns_for_file(filepath, args):
    """Return (smiles_col, id_col) for a given file, checking --column-map first."""
    basename = os.path.basename(filepath)
    if args.column_map:
        try:
            cmap = json.loads(args.column_map)
        except json.JSONDecodeError:
            cmap = {}
        if basename in cmap:
            entry = cmap[basename]
            return entry.get("smiles", args.smiles_col), entry.get("id", args.id_col)
    return args.smiles_col, args.id_col


# ── Molecule iterators ───────────────────────────────────────────────────────

def iter_csv(filepath, smiles_col, id_col):
    """Yield (mol_id, smiles, rdkit_mol_or_None, error_or_None) from a CSV."""
    with open(filepath, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, start=1):
            mol_id = row.get(id_col, f"row_{idx}")
            smiles = row.get(smiles_col, "").strip()
            if not smiles:
                yield mol_id, "", None, "Empty SMILES"
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                yield mol_id, smiles, None, "RDKit could not parse SMILES"
            else:
                yield mol_id, smiles, mol, None


def iter_sdf(filepath):
    """Yield (mol_id, smiles, rdkit_mol_or_None, error_or_None) from an SDF."""
    with open(filepath, "rb") as f:
        supplier = Chem.ForwardSDMolSupplier(f)
        for idx, mol in enumerate(supplier):
            if mol is None:
                yield f"mol_{idx}", "", None, f"Failed to parse molecule at index {idx}"
                continue
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
            mol_id = name.strip() or f"mol_{idx}"
            try:
                smiles = Chem.MolToSmiles(mol)
            except Exception as e:
                yield mol_id, "", None, f"Could not generate SMILES: {e}"
                continue
            yield mol_id, smiles, mol, None


def collect_input_files(input_path):
    """Return a list of (filepath, format) tuples. Supports file or folder."""
    if os.path.isfile(input_path):
        ext = os.path.splitext(input_path)[1].lower()
        if ext in (".csv", ".tsv"):
            return [(input_path, "csv")]
        elif ext in (".sdf", ".sd", ".mol"):
            return [(input_path, "sdf")]
        else:
            sys.exit(f"Error: unsupported file type '{ext}'. Use .csv or .sdf")
    elif os.path.isdir(input_path):
        files = []
        for pattern, fmt in [("*.csv", "csv"), ("*.tsv", "csv"), ("*.sdf", "sdf"), ("*.sd", "sdf")]:
            for fp in sorted(glob.glob(os.path.join(input_path, pattern))):
                files.append((fp, fmt))
        if not files:
            sys.exit(f"Error: no CSV or SDF files found in '{input_path}'")
        return files
    else:
        sys.exit(f"Error: '{input_path}' is not a file or directory")


# ── Standardization + writing ────────────────────────────────────────────────

def process_molecule(mol_id, smiles, mol, parse_error, source_file, options, writer, fail_writer):
    """Standardize one molecule and write to the appropriate output(s).
    Returns (is_success: bool)."""
    basename = os.path.basename(source_file)

    if parse_error:
        writer.writerow({
            ID_COL: mol_id, "source_file": basename,
            "original_smiles": smiles, "standardized_smiles": "",
            "success": False, "error": parse_error,
        })
        fail_writer.writerow({
            ID_COL: mol_id, "source_file": basename,
            "original_smiles": smiles, "error": parse_error, "checker_issues": "",
        })
        return False

    result = standardize_molecule(mol, options)

    checker_str = (
        "; ".join(f"[{score}] {msg}" for score, msg in result.checker_issues)
        if result.checker_issues else ""
    )

    writer.writerow({
        ID_COL: mol_id,
        "source_file": basename,
        "original_smiles": result.original_smiles,
        "standardized_smiles": result.standardized_smiles or "",
        "success": result.success,
        "error": result.error_message or "",
        "checker_issues": checker_str,
        "excluded_fragments": (
            "; ".join(result.excluded_fragments) if result.excluded_fragments else ""
        ),
        "stereo_warning": (
            result.stereo_comparison.warning
            if result.stereo_comparison and result.stereo_comparison.warning else ""
        ),
        "mass_change_pct": f"{result.mass_change_percent:.2f}" if result.success else "",
        "formula_original": (
            result.structure_comparison.original_formula
            if result.structure_comparison else ""
        ),
        "formula_standardized": (
            result.structure_comparison.standardized_formula
            if result.structure_comparison else ""
        ),
        "steps_summary": " -> ".join(
            f"{s.step_name}({'applied' if s.applied else 'skipped'})"
            for s in result.steps_applied
        ),
    })

    if not result.success:
        fail_writer.writerow({
            ID_COL: mol_id, "source_file": basename,
            "original_smiles": result.original_smiles,
            "error": result.error_message or "Unknown error",
            "checker_issues": checker_str,
        })

    return result.success



def main():
    args = parse_args()

    input_files = collect_input_files(args.input)
    is_folder = os.path.isdir(args.input)

    print(f"Found {len(input_files)} input file(s):")
    for fp, fmt in input_files:
        print(f"  [{fmt.upper()}] {fp}")
    print()

    # --output is treated as a directory when input is a folder with multiple files,
    # or as a single file path when input is a single file.
    if is_folder and len(input_files) > 1:
        out_dir = args.output
    else:
        out_dir = os.path.dirname(args.output) or "."
    os.makedirs(out_dir, exist_ok=True)

    options = StandardizationOptions(
        include_tautomer=args.include_tautomer,
        preserve_stereo=True,
        return_excluded_fragments=True,
    )

    t0 = time.time()
    grand_total = 0
    grand_success = 0
    grand_fail = 0

    # Build list of (input_path, fmt, output_path, failed_path)
    file_plan = []
    if is_folder and len(input_files) > 1:
        for filepath, fmt in input_files:
            stem = os.path.splitext(os.path.basename(filepath))[0]
            out_path = os.path.join(out_dir, f"{stem}_standardized.csv")
            fail_path = os.path.join(out_dir, f"{stem}_failed.csv")
            file_plan.append((filepath, fmt, out_path, fail_path))
    else:
        # Single file — use --output directly
        fp, fmt = input_files[0]
        base, ext = os.path.splitext(args.output)
        fail_path = f"{base}_failed{ext}"
        file_plan.append((fp, fmt, args.output, fail_path))

    for filepath, fmt, out_path, fail_path in file_plan:
        basename = os.path.basename(filepath)
        print(f"Processing {basename}...")
        file_total = 0
        file_success = 0
        file_fail = 0

        with open(out_path, "w", newline="", encoding="utf-8") as fout, \
             open(fail_path, "w", newline="", encoding="utf-8") as ffail:

            writer = csv.DictWriter(fout, fieldnames=FIELDNAMES, extrasaction="ignore")
            writer.writeheader()
            fail_writer = csv.DictWriter(ffail, fieldnames=FAILED_FIELDNAMES, extrasaction="ignore")
            fail_writer.writeheader()

            if fmt == "csv":
                smiles_col, id_col = get_columns_for_file(filepath, args)
                mol_iter = iter_csv(filepath, smiles_col, id_col)
            else:
                mol_iter = iter_sdf(filepath)

            for mol_id, smiles, mol, err in mol_iter:
                file_total += 1
                ok = process_molecule(
                    mol_id, smiles, mol, err, filepath, options, writer, fail_writer,
                )
                if ok:
                    file_success += 1
                else:
                    file_fail += 1

                if file_total % 200 == 0:
                    print(f"  Processed {file_total} molecules...")

        grand_total += file_total
        grand_success += file_success
        grand_fail += file_fail

        print(f"  {basename}: {file_total} total, {file_success} ok, {file_fail} failed")
        print(f"    -> {out_path}")
        if file_fail > 0:
            print(f"    -> {fail_path}")
        else:
            # Remove empty failed file
            os.remove(fail_path)
        print()

    elapsed = time.time() - t0
    print(f"Done. {grand_total} molecules across {len(file_plan)} file(s) in {elapsed:.1f}s")
    print(f"  Success: {grand_success}  |  Failed: {grand_fail}")



if __name__ == "__main__":
    main()

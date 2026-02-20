#!/bin/bash
#
# Interactive SLURM script generator for ChemAudit standardization.
# Queries live SLURM node availability (same strategy as run_orca_multi.sh),
# lets you pick a partition and auto-selects the best node.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# ── Helper: prompt for a path with Tab completion ────────────────────────────
# Uses readline (-e) so Tab completes files/dirs in the terminal.
# Args: $1 = prompt text, $2 = default value
# Sets: REPLY_PATH
read_path() {
    local prompt="$1" default="$2"
    # Save current completion setting, switch to filename completion
    local old_comp
    old_comp=$(bind -v 2>/dev/null | grep 'complete-filenames' || true)
    bind 'set show-all-if-ambiguous on' 2>/dev/null || true

    # -e enables readline (Tab completion), -i pre-fills the default
    read -rep "$prompt" -i "$default" REPLY_PATH
    REPLY_PATH="${REPLY_PATH:-$default}"
}

echo "=== ChemAudit Standardization — SLURM Script Generator ==="
echo "(Tip: path prompts support Tab completion)"
echo ""

# ── Input (file or folder) ───────────────────────────────────────────────────
DEFAULT_INPUT="$PROJECT_ROOT/examples/260219_Antimetabolite SMILES_H37Rv_Master.csv"
echo "Input can be a single CSV/SDF file, or a folder containing .csv/.sdf files."
read_path "Input file or folder: " "$DEFAULT_INPUT"
INPUT_PATH="$REPLY_PATH"

if [[ -d "$INPUT_PATH" ]]; then
    FILE_COUNT=$(find "$INPUT_PATH" -maxdepth 1 -type f \( -iname '*.csv' -o -iname '*.sdf' -o -iname '*.tsv' -o -iname '*.sd' \) 2>/dev/null | wc -l)
    echo "Folder mode: found $FILE_COUNT CSV/SDF file(s) in $INPUT_PATH"
    IS_FOLDER=true
elif [[ -f "$INPUT_PATH" ]]; then
    echo "Single file mode: $INPUT_PATH"
    IS_FOLDER=false
else
    echo "Warning: '$INPUT_PATH' not found locally. Continuing anyway (may exist on cluster)."
    IS_FOLDER=false
fi

if [[ "$IS_FOLDER" == true ]]; then
    DEFAULT_OUTPUT="$PROJECT_ROOT/results"
    echo "Each input file will get its own <name>_standardized.csv and <name>_failed.csv."
    read_path "Output directory: " "$DEFAULT_OUTPUT"
else
    DEFAULT_OUTPUT="$PROJECT_ROOT/results/standardized_output.csv"
    read_path "Output CSV path: " "$DEFAULT_OUTPUT"
fi
OUTPUT_CSV="$REPLY_PATH"

# ── Per-file column selection (CSV only) ─────────────────────────────────────
# Detect headers in each CSV and let the user pick SMILES and ID columns.
# SDF files use molecule name (_Name) and structure directly — no columns needed.

# Helper: present numbered column list and let user choose
# Args: $1 = prompt, $2 = columns (newline-separated), $3 = auto-detect hint (optional)
# Sets: CHOSEN_COL
pick_column() {
    local prompt="$1" col_str="$2" hint="${3:-}"
    local -a cols
    IFS=$'\n' read -rd '' -a cols <<< "$col_str" || true

    # Try to auto-detect a sensible default
    local default_idx=""
    if [[ -n "$hint" ]]; then
        for i in "${!cols[@]}"; do
            # Case-insensitive match
            if [[ "${cols[$i],,}" == *"${hint,,}"* ]]; then
                default_idx=$((i + 1))
                break
            fi
        done
    fi

    for i in "${!cols[@]}"; do
        local marker=""
        if [[ -n "$default_idx" ]] && [[ $((i + 1)) -eq $default_idx ]]; then
            marker=" <--"
        fi
        echo "    $((i+1))) ${cols[$i]}${marker}"
    done

    local sel_prompt="$prompt"
    if [[ -n "$default_idx" ]]; then
        sel_prompt+=" [$default_idx]"
    fi
    sel_prompt+=": "

    read -rp "$sel_prompt" SEL
    SEL="${SEL:-$default_idx}"

    if ! [[ "$SEL" =~ ^[0-9]+$ ]] || [ "$SEL" -lt 1 ] || [ "$SEL" -gt ${#cols[@]} ]; then
        echo "  Invalid selection, using column 1."
        SEL=1
    fi
    CHOSEN_COL="${cols[$((SEL-1))]}"
}

# Collect CSV files to inspect
declare -a CSV_FILES=()
if [[ "$IS_FOLDER" == true ]] && [[ -d "$INPUT_PATH" ]]; then
    while IFS= read -r -d '' f; do
        CSV_FILES+=("$f")
    done < <(find "$INPUT_PATH" -maxdepth 1 -type f \( -iname '*.csv' -o -iname '*.tsv' \) -print0 2>/dev/null | sort -z)
elif [[ -f "$INPUT_PATH" ]]; then
    ext="${INPUT_PATH##*.}"
    if [[ "${ext,,}" == "csv" || "${ext,,}" == "tsv" ]]; then
        CSV_FILES+=("$INPUT_PATH")
    fi
fi

# Build column map JSON: {"filename": {"smiles": "COL", "id": "COL"}, ...}
declare -A COL_MAP_SMILES
declare -A COL_MAP_ID

if [[ ${#CSV_FILES[@]} -gt 0 ]]; then
    echo ""
    echo "--- Column selection ---"
    for csv_file in "${CSV_FILES[@]}"; do
        fname="$(basename "$csv_file")"
        # Read the header line
        header_line=$(head -n1 "$csv_file" 2>/dev/null || true)
        if [[ -z "$header_line" ]]; then
            echo "  Warning: could not read headers from $fname, using defaults."
            COL_MAP_SMILES["$fname"]="SMILES"
            COL_MAP_ID["$fname"]="RegistrationNo"
            continue
        fi

        # Parse columns (handle comma and tab delimiters)
        local_delim=","
        if [[ "$header_line" == *$'\t'* ]]; then
            local_delim=$'\t'
        fi
        IFS="$local_delim" read -ra HEADERS <<< "$header_line"
        # Trim whitespace and quotes from each header
        CLEAN_HEADERS=()
        for h in "${HEADERS[@]}"; do
            h="$(echo "$h" | sed 's/^[[:space:]"]*//;s/[[:space:]"]*$//')"
            CLEAN_HEADERS+=("$h")
        done

        COL_LIST=$(printf '%s\n' "${CLEAN_HEADERS[@]}")

        echo ""
        echo "  File: $fname"
        echo "  Columns found:"

        pick_column "  Select SMILES column" "$COL_LIST" "smiles"
        COL_MAP_SMILES["$fname"]="$CHOSEN_COL"

        pick_column "  Select ID column" "$COL_LIST" "id"
        COL_MAP_ID["$fname"]="$CHOSEN_COL"

        echo "  -> SMILES='${COL_MAP_SMILES[$fname]}', ID='${COL_MAP_ID[$fname]}'"
    done
fi

# Build the --column-map JSON string
COLUMN_MAP_JSON="{"
FIRST_ENTRY=true
for fname in "${!COL_MAP_SMILES[@]}"; do
    if [[ "$FIRST_ENTRY" == true ]]; then
        FIRST_ENTRY=false
    else
        COLUMN_MAP_JSON+=","
    fi
    smiles_val="${COL_MAP_SMILES[$fname]}"
    id_val="${COL_MAP_ID[$fname]}"
    # Escape any double quotes in column names
    smiles_val="${smiles_val//\"/\\\"}"
    id_val="${id_val//\"/\\\"}"
    COLUMN_MAP_JSON+="\"$fname\":{\"smiles\":\"$smiles_val\",\"id\":\"$id_val\"}"
done
COLUMN_MAP_JSON+="}"

# If no CSV files, use empty map (SDF-only)
if [[ ${#CSV_FILES[@]} -eq 0 ]]; then
    COLUMN_MAP_JSON="{}"
fi

# ── Tautomer canonicalization ────────────────────────────────────────────────
read -rp "Enable tautomer canonicalization? (may lose E/Z stereo) [y/N]: " TAUTOMER
TAUTOMER_FLAG=""
if [[ "${TAUTOMER,,}" == "y" || "${TAUTOMER,,}" == "yes" ]]; then
    TAUTOMER_FLAG="--include-tautomer"
fi

# ── Pixi module ─────────────────────────────────────────────────────────────
echo ""
read -rp "Pixi module to load [pixi/0.56.0]: " PIXI_MODULE
PIXI_MODULE="${PIXI_MODULE:-pixi/0.56.0}"

# ── CPU request ──────────────────────────────────────────────────────────────
echo ""
echo "--- SLURM resource selection ---"
echo ""

read -rp "CPUs needed [4]: " CPUS_PER_NODE
CPUS_PER_NODE="${CPUS_PER_NODE:-4}"
if ! [[ "$CPUS_PER_NODE" =~ ^[0-9]+$ ]] || [ "$CPUS_PER_NODE" -lt 1 ]; then
    echo "Error: please enter a positive integer."
    exit 1
fi

echo ""
echo "Searching for nodes with at least $CPUS_PER_NODE free CPUs..."
echo ""

# ── Query SLURM for usable nodes ────────────────────────────────────────────
declare -A PART_NODES  # partition -> list of "node:total:alloc:free:free_mem:state"

while read -r part node total_cpus state; do
    case "$state" in
        idle|mix|idle*|mix*) ;;
        *) continue ;;
    esac
    NODE_INFO=$(scontrol show node "$node" 2>/dev/null)
    alloc_cpus=$(echo "$NODE_INFO" | grep -oP 'CPUAlloc=\K[0-9]+' || echo "0")
    free_cpus=$(( total_cpus - alloc_cpus ))
    real_mem=$(echo "$NODE_INFO" | grep -oP 'RealMemory=\K[0-9]+' || echo "0")
    alloc_mem=$(echo "$NODE_INFO" | grep -oP 'AllocMem=\K[0-9]+' || echo "0")
    free_mem=$(( real_mem - alloc_mem ))
    if [ "$free_cpus" -ge "$CPUS_PER_NODE" ]; then
        clean_part="${part%\*}"
        PART_NODES["$clean_part"]+="${node}:${total_cpus}:${alloc_cpus}:${free_cpus}:${free_mem}:${state} "
    fi
done < <(sinfo -h -N -o "%P %N %c %t" 2>/dev/null)

# ── Find partitions with at least 1 suitable node ───────────────────────────
declare -a VALID_PARTS=()
for part in "${!PART_NODES[@]}"; do
    count=$(echo "${PART_NODES[$part]}" | wc -w)
    if [ "$count" -ge 1 ]; then
        VALID_PARTS+=("$part")
    fi
done

if [ ${#VALID_PARTS[@]} -eq 0 ]; then
    echo "Error: no partition has a node with at least $CPUS_PER_NODE free CPUs."
    echo "(Draining, drained, down, and reserved nodes are excluded.)"
    exit 1
fi

# ── Let user pick a partition ────────────────────────────────────────────────
echo "Partitions with suitable nodes:"
echo ""
for pi in "${!VALID_PARTS[@]}"; do
    part="${VALID_PARTS[$pi]}"
    count=$(echo "${PART_NODES[$part]}" | wc -w)
    echo "  $((pi+1))) $part ($count node(s) available)"
done
echo ""

if [ ${#VALID_PARTS[@]} -eq 1 ]; then
    CHOSEN_PARTITION="${VALID_PARTS[0]}"
    echo "Auto-selected partition: $CHOSEN_PARTITION"
else
    read -rp "Select a partition [1-${#VALID_PARTS[@]}]: " PSEL
    if ! [[ "$PSEL" =~ ^[0-9]+$ ]] || [ "$PSEL" -lt 1 ] || [ "$PSEL" -gt ${#VALID_PARTS[@]} ]; then
        echo "Error: invalid selection."
        exit 1
    fi
    CHOSEN_PARTITION="${VALID_PARTS[$((PSEL-1))]}"
fi

echo ""

# ── Auto-select best node: prefer idle, then most free memory ────────────────
TMPFILE=$(mktemp)
for entry in ${PART_NODES[$CHOSEN_PARTITION]}; do
    IFS=':' read -r node total alloc free free_mem state <<< "$entry"
    case "$state" in
        idle|idle*) prio=0 ;;
        *)          prio=1 ;;
    esac
    echo "${prio}:${free_mem}:${node}:${free}:${state}" >> "$TMPFILE"
done

sort -t: -k1,1n -k2,2rn "$TMPFILE" -o "$TMPFILE"

echo "Best available nodes on '$CHOSEN_PARTITION':"
echo ""
printf "  %-4s %-12s %-10s %-10s %s\n" "#" "NODE" "FREE_CPUs" "FREE_MEM" "STATE"
printf "  %-4s %-12s %-10s %-10s %s\n" "---" "----------" "---------" "--------" "-----"

IDX=0
declare -a NODE_LIST=()
declare -a NODE_FREEMEM=()
while IFS=: read -r prio fmem node fcpus state; do
    IDX=$(( IDX + 1 ))
    NODE_LIST+=("$node")
    NODE_FREEMEM+=("$fmem")
    fmem_gb=$(( fmem / 1024 ))
    printf "  %-4s %-12s %-10s %-10s %s\n" "$IDX)" "$node" "$fcpus" "${fmem_gb}GB" "$state"
done < "$TMPFILE"

rm -f "$TMPFILE"

echo ""

if [ ${#NODE_LIST[@]} -eq 1 ]; then
    CHOSEN_NODE="${NODE_LIST[0]}"
    CHOSEN_FREE_MEM="${NODE_FREEMEM[0]}"
    echo "Auto-selected node: $CHOSEN_NODE"
else
    echo "  0) Any node (let SLURM decide)"
    echo ""
    read -rp "Select a node [0-${#NODE_LIST[@]}, default=1]: " NSEL
    NSEL="${NSEL:-1}"

    if [[ "$NSEL" == "0" ]]; then
        CHOSEN_NODE=""
        CHOSEN_FREE_MEM=""
    elif [[ "$NSEL" =~ ^[0-9]+$ ]] && [ "$NSEL" -ge 1 ] && [ "$NSEL" -le ${#NODE_LIST[@]} ]; then
        CHOSEN_NODE="${NODE_LIST[$((NSEL-1))]}"
        CHOSEN_FREE_MEM="${NODE_FREEMEM[$((NSEL-1))]}"
    else
        echo "Error: invalid selection."
        exit 1
    fi
fi

# ── Calculate memory ─────────────────────────────────────────────────────────
if [[ -n "$CHOSEN_FREE_MEM" ]]; then
    USABLE_MEM=$(( CHOSEN_FREE_MEM * 95 / 100 ))
    MEM_STR="${USABLE_MEM}M"
    echo "Free memory on node: ${CHOSEN_FREE_MEM}MB → requesting ${USABLE_MEM}MB (95%)"
else
    read -rp "Memory to request (e.g. 8G) [8G]: " MEM_STR
    MEM_STR="${MEM_STR:-8G}"
fi

# ── Remaining SLURM options ──────────────────────────────────────────────────
echo ""
read -rp "Wall time (HH:MM:SS or D-HH:MM:SS) [01:00:00]: " WALL_TIME
WALL_TIME="${WALL_TIME:-01:00:00}"

read -rp "Job name [chemaudit-std]: " JOB_NAME
JOB_NAME="${JOB_NAME:-chemaudit-std}"

read -rp "Email for notifications (leave blank to skip): " EMAIL

# ── Output script path ──────────────────────────────────────────────────────
echo ""
DEFAULT_SCRIPT="$SCRIPT_DIR/submit_standardization.sh"
read_path "Output SLURM script path: " "$DEFAULT_SCRIPT"
SLURM_SCRIPT="$REPLY_PATH"

# ── Build optional SBATCH lines ─────────────────────────────────────────────
NODELIST_LINE=""
if [[ -n "$CHOSEN_NODE" ]]; then
    NODELIST_LINE="#SBATCH --nodelist=$CHOSEN_NODE"
fi

EMAIL_LINES=""
if [[ -n "$EMAIL" ]]; then
    EMAIL_LINES="#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=$EMAIL"
fi

# Compute the mkdir target for the generated script
if [[ "$IS_FOLDER" == true ]]; then
    MKDIR_TARGET="$OUTPUT_CSV"
else
    MKDIR_TARGET="$(dirname "$OUTPUT_CSV")"
fi

# ── Write the SLURM script ──────────────────────────────────────────────────
cat > "$SLURM_SCRIPT" <<SLURM_EOF
#!/bin/bash
#SBATCH --job-name=$JOB_NAME
#SBATCH --output=results/slurm-%j.out
#SBATCH --error=results/slurm-%j.err
#SBATCH --time=$WALL_TIME
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$CPUS_PER_NODE
#SBATCH --mem=$MEM_STR
#SBATCH --partition=$CHOSEN_PARTITION
${NODELIST_LINE}
${EMAIL_LINES}

# ── Load pixi ────────────────────────────────────────────────────────────────
module load $PIXI_MODULE

# ── Paths (hardcoded at generation time — SLURM copies scripts to /tmp) ──────
SCRIPT_DIR="$SCRIPT_DIR"

# ── Install dependencies on first run ────────────────────────────────────────
echo "Installing pixi environment..."
pixi install --manifest-path "\$SCRIPT_DIR/pixi.toml"

# ── Create output directory ───────────────────────────────────────────────────
mkdir -p "$MKDIR_TARGET"

# ── Run the standardization pipeline ─────────────────────────────────────────
echo "Starting standardization pipeline..."
echo "  Input:  $INPUT_PATH"
echo "  Output: $OUTPUT_CSV"
echo ""

pixi run --manifest-path "\$SCRIPT_DIR/pixi.toml" \\
    python "\$SCRIPT_DIR/run_standardization.py" \\
    --input "$INPUT_PATH" \\
    --output "$OUTPUT_CSV" \\
    --column-map '$COLUMN_MAP_JSON' $TAUTOMER_FLAG

echo ""
echo "Job finished at \$(date)"
SLURM_EOF

chmod +x "$SLURM_SCRIPT"

# ── Clean up blank lines from optional SBATCH directives ────────────────────
sed -i '/^$/{ N; /^\n$/d; }' "$SLURM_SCRIPT"

echo ""
echo "=== Generated: $SLURM_SCRIPT ==="
echo "---"
cat "$SLURM_SCRIPT"
echo "---"
echo ""

# ── Offer to submit ─────────────────────────────────────────────────────────
read -rp "Submit now? [y/N]: " CONFIRM
if [[ "$CONFIRM" =~ ^[Yy]$ ]]; then
    mkdir -p "$PROJECT_ROOT/results"
    sbatch "$SLURM_SCRIPT"
else
    echo "Not submitted. Run manually with: sbatch $SLURM_SCRIPT"
fi

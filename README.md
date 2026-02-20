<div align="center">

<img src="docs/assets/Banner_logo.png" alt="ChemAudit" width="600" />

### Chemical Structure Validation & Quality Assessment Platform

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18492676.svg)](https://doi.org/10.5281/zenodo.18492676)
[![Release](https://img.shields.io/github/v/release/Kohulan/ChemAudit)](https://github.com/Kohulan/ChemAudit/releases)
[![Tests](https://github.com/Kohulan/ChemAudit/actions/workflows/test.yml/badge.svg)](https://github.com/Kohulan/ChemAudit/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen.svg)](https://kohulan.github.io/ChemAudit/)
[![Contributors](https://img.shields.io/github/contributors/Kohulan/ChemAudit)](https://github.com/Kohulan/ChemAudit/graphs/contributors)
[![Issues](https://img.shields.io/github/issues/Kohulan/ChemAudit)](https://github.com/Kohulan/ChemAudit/issues)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![React 18](https://img.shields.io/badge/react-18-61DAFB.svg?logo=react)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-009688.svg?logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![RDKit](https://img.shields.io/badge/RDKit-3776AB.svg)](https://www.rdkit.org/)
[![Docker](https://img.shields.io/badge/docker-ready-2496ED.svg?logo=docker&logoColor=white)](https://www.docker.com/)

<br />

**Validate • Standardize • Score • Analyze**

*A powerful web platform for cheminformatics workflows, drug discovery, and ML dataset curation*

[Features](#-features) •
[Quick Start](#-quick-start) •
[Documentation](https://www.kohulanr.com/ChemAudit/) •
[API](#-api-reference) •
[Contributing](#-contributing)

<br />

<img src="docs/assets/screenshot-dashboard.png" alt="ChemAudit Dashboard" width="800" />

</div>

---

## ✨ Features

<table>
<tr>
<td width="50%">

### 🔬 Structure Validation
Comprehensive chemical structure analysis with 15+ validation checks

- Valence & connectivity errors
- Aromaticity & kekulization
- Stereochemistry validation
- Ring system analysis
- Atom & bond type verification

</td>
<td width="50%">

### ⚠️ Structural Alerts
Screen compounds against known problematic substructures

- **PAINS** - Pan-Assay Interference detection
- **BRENK** - Unwanted group filters
- **Custom Catalogs** - Define your own rules
- Severity scoring & recommendations

</td>
</tr>
<tr>
<td width="50%">

### 📊 ML-Readiness Scoring
Evaluate compound suitability for machine learning

- Descriptor calculability assessment
- Fingerprint generation validation
- Molecular complexity scoring
- Dataset quality metrics

</td>
<td width="50%">

### 🌿 Natural Product Analysis
Assess natural product-likeness characteristics

- NP-Likeness scoring algorithm
- Scaffold extraction & classification
- Structural similarity to known NPs
- Biosynthetic pathway hints

</td>
</tr>
<tr>
<td width="50%">

### 🧹 Standardization Pipeline
ChEMBL-compatible molecular standardization

- Salt stripping & neutralization
- Tautomer canonicalization
- Stereochemistry normalization
- Configurable pipeline steps

</td>
<td width="50%">

### 🗄️ Database Integrations
Cross-reference against major chemical databases

- **PubChem** - Compound properties & synonyms
- **ChEMBL** - Bioactivity & target data
- **COCONUT** - Natural products database

</td>
</tr>
</table>

---

## 🚀 Quick Start

### Using Docker (Recommended)

```bash
# Clone the repository
git clone https://github.com/Kohulan/ChemAudit.git
cd chemaudit

# Start all services (development)
docker-compose up -d

# View logs
docker-compose logs -f
```

🌐 **Access Points:**
| Service | URL |
|---------|-----|
| Web UI | http://localhost:3002 |
| API Docs | http://localhost:8001/docs |
| Metrics | http://localhost:9090 |

### Production Deployment

Use the interactive deploy script to select a deployment profile:

```bash
# Interactive mode - shows profile menu
./deploy.sh

# Or specify profile directly
./deploy.sh medium
```

**Available Profiles:**

| Profile | Max Molecules | Max File Size | Workers | Use Case |
|---------|---------------|---------------|---------|----------|
| `small` | 1,000 | 100 MB | 2 | Development |
| `medium` | 10,000 | 500 MB | 4 | Standard production |
| `large` | 50,000 | 500 MB | 8 | High-throughput |
| `xl` | 100,000 | 1 GB | 12 | Enterprise |
| `coconut` | 1,000,000 | 1 GB | 16 | Full COCONUT DB |

See [Deployment Guide](docs/DEPLOYMENT.md) for detailed configuration.

### Manual Installation

<details>
<summary><b>Backend Setup</b></summary>

```bash
cd backend
poetry install
poetry run uvicorn app.main:app --reload
```

</details>

<details>
<summary><b>Frontend Setup</b></summary>

```bash
cd frontend
npm install
npm run dev
```

</details>

---

## 📸 Screenshots

<div align="center">
<table>
<tr>
<td align="center"><b>Single Molecule Validation</b></td>
<td align="center"><b>Batch Processing</b></td>
</tr>
<tr>
<td><img src="docs/assets/screenshot-single.png" width="400" /></td>
<td><img src="docs/assets/screenshot-batch.png" width="400" /></td>
</tr>
<tr>
<td align="center"><b>Scoring Dashboard</b></td>
<td align="center"><b>Database Lookup</b></td>
</tr>
<tr>
<td><img src="docs/assets/screenshot-scoring.png" width="400" /></td>
<td><img src="docs/assets/screenshot-integrations.png" width="400" /></td>
</tr>
</table>
</div>

---

## 📦 Batch Processing

Process large datasets with ease:

| Feature | Specification |
|---------|---------------|
| **Max File Size** | Up to 1 GB (profile-dependent) |
| **Max Molecules** | Up to 1,000,000 per batch (profile-dependent) |
| **Supported Formats** | SDF, CSV |
| **Progress Tracking** | Real-time WebSocket updates |
| **Export Formats** | CSV, JSON, Excel, SDF, PDF Report |

> **Note:** Batch limits are dynamically configured based on your deployment profile.
> The frontend automatically displays the correct limits for your deployment.

```python
# Python client example
from chemaudit import ChemAuditClient

client = ChemAuditClient("http://localhost:8000")

# Upload and process
job = client.upload_batch("molecules.sdf")

# Monitor progress
for update in client.stream_progress(job.job_id):
    print(f"Progress: {update.progress}%")

# Get results
results = client.get_results(job.job_id)
```

---

## 🛠️ Tech Stack

<div align="center">

| Layer | Technologies |
|-------|-------------|
| **Frontend** | ![React](https://img.shields.io/badge/React-61DAFB?logo=react&logoColor=black) ![TypeScript](https://img.shields.io/badge/TypeScript-3178C6?logo=typescript&logoColor=white) ![Vite](https://img.shields.io/badge/Vite-646CFF?logo=vite&logoColor=white) ![Tailwind](https://img.shields.io/badge/Tailwind-06B6D4?logo=tailwindcss&logoColor=white) |
| **Backend** | ![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=white) ![FastAPI](https://img.shields.io/badge/FastAPI-009688?logo=fastapi&logoColor=white) ![RDKit](https://img.shields.io/badge/RDKit-3776AB) ![Celery](https://img.shields.io/badge/Celery-37814A?logo=celery&logoColor=white) |
| **Database** | ![PostgreSQL](https://img.shields.io/badge/PostgreSQL-4169E1?logo=postgresql&logoColor=white) ![Redis](https://img.shields.io/badge/Redis-DC382D?logo=redis&logoColor=white) |
| **Infrastructure** | ![Docker](https://img.shields.io/badge/Docker-2496ED?logo=docker&logoColor=white) ![Nginx](https://img.shields.io/badge/Nginx-009639?logo=nginx&logoColor=white) |
| **Monitoring** | ![Prometheus](https://img.shields.io/badge/Prometheus-E6522C?logo=prometheus&logoColor=white) ![Grafana](https://img.shields.io/badge/Grafana-F46800?logo=grafana&logoColor=white) |

</div>

---

## 📖 Documentation

| Document | Description |
|----------|-------------|
| [Getting Started](docs/GETTING_STARTED.md) | Installation and first steps |
| [User Guide](docs/USER_GUIDE.md) | Complete usage instructions |
| [API Reference](docs/API_REFERENCE.md) | Full REST API documentation |
| [Deployment](docs/DEPLOYMENT.md) | Production deployment guide |
| [Troubleshooting](docs/TROUBLESHOOTING.md) | Common issues & solutions |

**Interactive API Docs:** http://localhost:8000/docs

---

## 🔌 API Reference

### Validate a Molecule

```bash
curl -X POST http://localhost:8000/api/v1/validate \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CC(=O)Oc1ccccc1C(=O)O", "format": "smiles"}'
```

<details>
<summary><b>Response</b></summary>

```json
{
  "valid": true,
  "validation_score": 92,
  "checks": [
    {"name": "valence", "passed": true, "severity": "critical"},
    {"name": "aromaticity", "passed": true, "severity": "warning"},
    {"name": "stereo", "passed": true, "severity": "info"}
  ],
  "standardized_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
  "inchi_key": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
}
```

</details>

### Screen for Alerts

```bash
curl -X POST http://localhost:8000/api/v1/alerts \
  -H "Content-Type: application/json" \
  -d '{"molecule": "c1ccc2c(c1)nc(n2)Sc3nnnn3C", "catalogs": ["PAINS", "BRENK"]}'
```

### Get ML-Readiness Score

```bash
curl -X POST http://localhost:8000/api/v1/score \
  -H "Content-Type: application/json" \
  -d '{"molecule": "CCO", "include": ["ml_readiness", "np_likeness"]}'
```

---

## 🏗️ Project Structure

```
chemaudit/
├── backend/
│   ├── app/
│   │   ├── api/routes/       # API endpoints
│   │   ├── core/             # Config, security, metrics
│   │   ├── schemas/          # Pydantic models
│   │   ├── services/         # Business logic
│   │   │   ├── validation/   # Validation checks
│   │   │   ├── scoring/      # ML-readiness, NP-likeness
│   │   │   ├── alerts/       # Structural alert screening
│   │   │   └── integrations/ # External DB connectors
│   │   └── templates/        # Report templates
│   └── tests/
├── frontend/
│   ├── src/
│   │   ├── components/       # React components
│   │   ├── pages/            # Route pages
│   │   ├── services/         # API client
│   │   ├── hooks/            # Custom hooks
│   │   └── types/            # TypeScript definitions
│   └── tests/
├── scripts/
│   ├── generate_slurm.sh    # Interactive SLURM generator
│   ├── run_standardization.py # CLI standardization pipeline
│   └── pixi.toml            # Pixi dependency manifest
├── nginx/                    # Reverse proxy config
├── docs/                     # Documentation
└── docker-compose.yml
```

---

## 🧪 Running Tests

```bash
# Backend tests
cd backend
poetry run pytest

# Frontend tests
cd frontend
npm test

# With coverage
poetry run pytest --cov=app
npm run test:coverage
```

---

## 🧹 CLI Standardization (HPC / Local)

Run the ChEMBL standardization pipeline directly on CSV or SDF files without the web server. Designed for HPC clusters with SLURM but works locally too.

### Files

| File | Purpose |
|------|---------|
| `scripts/run_standardization.py` | Standalone standardization script |
| `scripts/generate_slurm.sh` | Interactive SLURM job generator |
| `scripts/pixi.toml` | Conda/PyPI dependency manifest |

### Local Usage

```bash
# Install dependencies with pixi (https://pixi.sh)
pixi install --manifest-path scripts/pixi.toml

# Single CSV file
pixi run --manifest-path scripts/pixi.toml \
    python scripts/run_standardization.py \
    --input data/molecules.csv \
    --output results/standardized.csv

# Single SDF file
pixi run --manifest-path scripts/pixi.toml \
    python scripts/run_standardization.py \
    --input data/molecules.sdf \
    --output results/standardized.csv

# Folder of CSV/SDF files (each gets its own output)
pixi run --manifest-path scripts/pixi.toml \
    python scripts/run_standardization.py \
    --input data/molecules/ \
    --output results/
```

When given a folder, each input file produces its own output pair:
- `mols_A.csv` → `results/mols_A_standardized.csv` + `results/mols_A_failed.csv`
- `mols_B.sdf` → `results/mols_B_standardized.csv` + `results/mols_B_failed.csv`

Failed molecule files are only created when there are actual failures.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | *(required)* | CSV/SDF file or folder |
| `--output` | *(required)* | Output file (single) or directory (folder) |
| `--smiles-col` | `SMILES` | Default SMILES column for CSVs |
| `--id-col` | `RegistrationNo` | Default ID column for CSVs |
| `--column-map` | — | Per-file column mapping as JSON (see below) |
| `--include-tautomer` | off | Enable tautomer canonicalization (may lose E/Z stereo) |

#### Per-file column mapping

When processing multiple CSVs with different column names, pass a JSON map:

```bash
--column-map '{"file_A.csv": {"smiles": "CANONICAL_SMILES", "id": "COMPOUND_ID"}, "file_B.csv": {"smiles": "smi", "id": "name"}}'
```

### HPC / SLURM Usage

The interactive generator handles everything — dependency setup, partition/node selection, and column mapping:

```bash
# Load pixi on the cluster
module load pixi/0.56.0

# Run the interactive generator
bash scripts/generate_slurm.sh
```

The generator will:
1. Ask for input file/folder (with Tab completion)
2. Detect CSV column headers and let you pick SMILES/ID columns per file
3. Query SLURM for available nodes and free resources
4. Let you choose a partition and node (or let SLURM decide)
5. Auto-calculate memory from the selected node
6. Write a ready-to-submit SLURM script and offer to `sbatch` it

---

## 🤝 Contributing

Contributions are welcome! Please read our contributing guidelines before submitting a PR.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- [RDKit](https://www.rdkit.org/) - Cheminformatics toolkit
- [ChEMBL](https://www.ebi.ac.uk/chembl/) - Bioactivity database
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - Chemical database
- [COCONUT](https://coconut.naturalproducts.net/) - Natural products database

---

<div align="center">

**[⬆ Back to Top](#-chemaudit)**

Made with ☕ by [Kohulan.R](https://github.com/Kohulan) at [Friedrich Schiller University Jena](https://www.uni-jena.de/)

</div>

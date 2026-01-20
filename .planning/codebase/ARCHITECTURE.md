# Architecture

**Analysis Date:** 2026-01-20

## Pattern

**Type:** Layered Monolith (Backend) + SPA (Frontend)

**Description:** Clean separation between React SPA frontend and FastAPI backend with async processing via Celery workers. Designed for horizontal scaling.

## Layers

**Frontend (React SPA):**
- `components/` - UI components (shadcn/ui based)
- `pages/` - Route-level page components
- `hooks/` - React Query hooks for API state
- `services/` - API client layer (Axios)
- `types/` - TypeScript type definitions

**Backend (FastAPI):**
- `api/routes/` - HTTP endpoint handlers
- `services/` - Business logic layer
- `models/` - SQLAlchemy ORM models
- `schemas/` - Pydantic request/response schemas
- `tasks/` - Celery async task definitions
- `db/` - Database connection and repositories

**Core Services:**
- `services/validation/` - ValidationEngine + check implementations
- `services/standardization/` - Molecular standardization pipeline
- `services/parser/` - Molecule format parsing
- `services/export/` - Result export (CSV, Excel, PDF)

## Data Flow

```
User Input (SMILES/SDF)
    │
    ▼
React Frontend ─────► FastAPI `/api/v1/validate`
    │                      │
    │                      ▼
    │                MoleculeParser.parse()
    │                      │
    │                      ▼
    │                ValidationEngine.validate()
    │                      │
    │           ┌──────────┼──────────┐
    │           ▼          ▼          ▼
    │      BasicChecks  StereoChecks  AlertEngine
    │           │          │          │
    │           └──────────┼──────────┘
    │                      ▼
    │              Calculate Scores
    │                      │
    │                      ▼
    ◄────────────── ValidationResponse
```

**Batch Flow:**
```
File Upload
    │
    ▼
POST /api/v1/batch ──► Create Job Record (PostgreSQL)
    │                       │
    │                       ▼
    │                 Celery Task Queued (Redis)
    │                       │
    │                       ▼
    │                 Worker: Process Chunks
    │                       │
    │                       ▼
    ◄──── WebSocket Progress Updates
    │
    ▼
GET /api/v1/jobs/{id}/results
```

## Key Abstractions

**ValidationEngine** (`backend/app/services/validation/engine.py`):
- Central orchestrator for all validation checks
- Plugin architecture via `register_check()`
- Calculates overall and ML-readiness scores

**BaseCheck** (`backend/app/services/validation/checks/base.py`):
- Abstract base class for all validation checks
- Interface: `run(mol: Mol) -> CheckResult`
- Categories: basic, stereo, alerts, representation

**StandardizationPipeline** (`backend/app/services/standardization/pipeline.py`):
- Configurable pipeline of standardization steps
- Steps: SaltStripper, Normalizer, TautomerCanonicalizer, etc.
- Based on ChEMBL structure pipeline

**AlertEngine** (`backend/app/services/validation/checks/alerts.py`):
- SMARTS pattern matching for structural alerts
- Supports: PAINS, BRENK, NIH, ZINC, ChEMBL alert sets
- Uses RDKit FilterCatalog

## Entry Points

**Backend:**
- `backend/app/main.py` - FastAPI application factory
- `backend/app/tasks/celery_app.py` - Celery worker entry

**Frontend:**
- `frontend/src/main.tsx` - React app entry
- `frontend/src/App.tsx` - Router and providers

**API Endpoints:**
- `POST /api/v1/validate` - Single molecule validation
- `POST /api/v1/batch` - Batch file upload
- `GET /api/v1/jobs/{id}` - Job status
- `GET /api/v1/jobs/{id}/results` - Batch results
- `WS /ws/jobs/{id}` - Real-time progress

## Cross-Cutting Concerns

**Caching (Redis):**
- Validation results by InChIKey
- Alert pattern definitions
- Job progress state

**Error Handling:**
- Pydantic validation for all inputs
- Custom exceptions with HTTP status mapping
- Graceful degradation for non-critical failures

**Logging:**
- Structured JSON logging
- Request correlation IDs
- Configurable log levels

---

*Architecture analysis: 2026-01-20*

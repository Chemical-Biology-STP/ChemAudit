# Technology Stack

**Analysis Date:** 2026-01-20

## Languages

**Primary:**
- Python 3.11+ - Backend API, validation engine, cheminformatics processing
- TypeScript 5.1+ - Frontend React application with strict mode

**Secondary:**
- SQL - PostgreSQL database queries and schema
- YAML - Docker Compose, GitHub Actions CI/CD configuration
- JSON - API schemas, alert pattern definitions

## Runtime

**Environment:**
- Node.js (latest LTS) - Frontend build and development
- Python 3.11+ - Backend runtime with async support
- Docker - Containerized deployment

**Package Manager:**
- Poetry - Python dependency management (`backend/pyproject.toml`)
- npm - Frontend dependency management (`frontend/package.json`)
- Lockfile: Both `poetry.lock` and `package-lock.json` expected

## Frameworks

**Core:**
- FastAPI ^0.100.0 - Backend REST API framework with async support
- React 18.2+ - Frontend UI library with functional components and hooks
- Pydantic ^2.0.0 - Request/response validation and settings management
- Vite ^4.4.0 - Frontend build tool and dev server

**Testing:**
- pytest ^7.4.0 - Backend unit and integration testing
- pytest-asyncio ^0.21.0 - Async test support
- pytest-cov ^4.1.0 - Code coverage
- Vitest ^0.34.0 - Frontend unit testing
- @testing-library/react ^14.0.0 - React component testing
- httpx ^0.24.0 - Async HTTP client for API testing

**Build/Dev:**
- Docker + Docker Compose - Local development and production deployment
- GitHub Actions - CI/CD pipelines
- Black ^23.0.0 - Python code formatter
- isort ^5.12.0 - Python import sorter
- mypy ^1.4.0 - Python static type checking
- ESLint + Prettier - TypeScript/React linting and formatting

## Key Dependencies

**Critical (Backend Chemistry):**
- rdkit ^2023.9.0 - Core cheminformatics library for molecule parsing, validation, fingerprints
- molvs ^0.1.1 - Molecule standardization and validation (built on RDKit)
- chembl-structure-pipeline ^1.2.0 - ChEMBL standardization, salt stripping, parent extraction

**Critical (Backend Infrastructure):**
- uvicorn ^0.22.0 - ASGI server for FastAPI
- celery ^5.3.0 - Distributed task queue for batch processing
- redis ^4.6.0 - Caching, job queue broker, session management
- asyncpg ^0.28.0 - Async PostgreSQL driver
- sqlalchemy ^2.0.0 - ORM and database toolkit
- alembic ^1.11.0 - Database migrations

**Critical (Frontend):**
- @rdkit/rdkit ^2023.9.1 - RDKit.js WebAssembly for client-side molecule rendering
- @tanstack/react-query ^4.32.0 - Server state management, caching
- axios ^1.4.0 - HTTP client for API calls
- react-router-dom ^6.14.0 - Client-side routing

**UI Components:**
- tailwindcss ^3.3.0 - Utility-first CSS framework
- shadcn/ui - Component library (Radix UI + Tailwind)
- lucide-react ^0.263.0 - Icon library
- recharts ^2.7.0 - Charting library for score visualizations
- react-dropzone ^14.2.0 - File upload drag-and-drop

**Utilities:**
- pandas ^2.0.0 - Data manipulation for batch processing
- numpy ^1.24.0 - Numerical operations (RDKit dependency)
- aiofiles ^23.0.0 - Async file operations
- python-multipart ^0.0.6 - File upload handling
- python-jose ^3.3.0 - JWT token handling
- slowapi ^0.1.8 - Rate limiting
- pydantic-settings ^2.0.0 - Environment configuration
- clsx ^2.0.0 - Class name utilities
- tailwind-merge ^1.14.0 - Tailwind class merging

## Configuration

**Environment:**
- `.env` files for local development
- pydantic-settings for type-safe configuration loading
- Key environment variables:
  - `DATABASE_URL` - PostgreSQL connection string
  - `REDIS_URL` - Redis connection string
  - `CELERY_BROKER_URL` - Celery broker (Redis)
  - `CORS_ORIGINS` - Allowed frontend origins
  - `MAX_BATCH_SIZE` - Maximum molecules per batch (default: 100,000)
  - `MAX_FILE_SIZE_MB` - Upload limit (default: 100MB)
  - `VITE_API_URL` - Backend API URL for frontend

**Build:**
- `backend/pyproject.toml` - Python project configuration, dependencies
- `frontend/package.json` - Node.js project configuration, dependencies
- `frontend/vite.config.ts` - Vite build configuration
- `frontend/tsconfig.json` - TypeScript compiler options
- `frontend/tailwind.config.js` - Tailwind CSS configuration
- `frontend/postcss.config.js` - PostCSS configuration
- `docker-compose.yml` - Development services
- `docker-compose.prod.yml` - Production deployment
- `.github/workflows/ci.yml` - CI pipeline

## Platform Requirements

**Development:**
- Docker Desktop or Docker Engine + Compose
- Python 3.11+ with pip or Poetry
- Node.js LTS (18+) with npm
- Git

**Production:**
- Linux containers (Docker)
- PostgreSQL 15+
- Redis 7+
- Load balancer (nginx/HAProxy)
- CDN for static assets (optional)

**Deployment Targets:**
- Docker containers (Kubernetes or Docker Swarm compatible)
- Horizontal scaling: multiple backend API instances
- Celery workers: scalable batch processing

---

*Stack analysis: 2026-01-20*

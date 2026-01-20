# Testing

**Analysis Date:** 2026-01-20

## Framework

**Backend:**
- Framework: pytest ^7.4.0
- Async support: pytest-asyncio ^0.21.0
- Coverage: pytest-cov ^4.1.0
- HTTP testing: httpx ^0.24.0

**Frontend:**
- Framework: Vitest ^0.34.0
- Component testing: @testing-library/react ^14.0.0
- DOM utilities: @testing-library/jest-dom

## Test Structure

**Backend:**
```
backend/tests/
├── conftest.py              # Shared fixtures
├── test_validation/
│   ├── test_engine.py       # ValidationEngine tests
│   ├── test_basic_checks.py # Basic check tests
│   ├── test_stereo_checks.py
│   └── test_alerts.py
├── test_api/
│   ├── test_validation_routes.py
│   ├── test_batch_routes.py
│   └── test_jobs_routes.py
├── test_standardization/
│   └── test_pipeline.py
└── fixtures/
    └── molecules.py         # Test molecule fixtures
```

**Frontend:**
```
frontend/src/
├── components/
│   └── validation/
│       ├── ValidationResults.tsx
│       └── ValidationResults.test.tsx  # Co-located
├── hooks/
│   ├── useValidation.ts
│   └── useValidation.test.ts
└── __tests__/               # Integration tests
    └── pages/
```

## Running Tests

**Backend:**
```bash
cd backend

# Run all tests
pytest

# With coverage
pytest --cov=app --cov-report=html

# Run specific file
pytest tests/test_validation/test_basic_checks.py -v

# Run specific test
pytest tests/test_validation/test_basic_checks.py::test_valence_check -v

# Run with markers
pytest -m "not slow"
```

**Frontend:**
```bash
cd frontend

# Run all tests
npm test

# Watch mode
npm test -- --watch

# Coverage
npm test -- --coverage

# Run specific file
npm test -- ValidationResults.test.tsx
```

## Mocking Patterns

**Backend - RDKit molecules:**
```python
# tests/fixtures/molecules.py
import pytest
from rdkit import Chem

@pytest.fixture
def valid_mol():
    """Valid ethanol molecule."""
    return Chem.MolFromSmiles("CCO")

@pytest.fixture
def invalid_valence_mol():
    """Molecule with valence error (pentavalent carbon)."""
    return Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)

@pytest.fixture
def multi_fragment_mol():
    """Molecule with disconnected fragments."""
    return Chem.MolFromSmiles("CCO.CC")
```

**Backend - Database:**
```python
# tests/conftest.py
import pytest
from sqlalchemy.ext.asyncio import create_async_engine, AsyncSession

@pytest.fixture
async def db_session():
    """Create test database session."""
    engine = create_async_engine("sqlite+aiosqlite:///:memory:")
    async with AsyncSession(engine) as session:
        yield session

@pytest.fixture
def mock_redis(mocker):
    """Mock Redis client."""
    return mocker.patch("app.core.cache.redis_client")
```

**Frontend - API calls:**
```typescript
// mocks/api.ts
import { vi } from 'vitest';

export const mockValidationApi = {
  validate: vi.fn(),
  getChecks: vi.fn(),
};

// In test file
vi.mock('@/services/api', () => ({
  validationApi: mockValidationApi,
}));

beforeEach(() => {
  mockValidationApi.validate.mockResolvedValue({
    overall_score: 85,
    issues: [],
    // ...
  });
});
```

**Frontend - React Query:**
```typescript
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render } from '@testing-library/react';

const queryClient = new QueryClient({
  defaultOptions: {
    queries: { retry: false },
  },
});

function renderWithProviders(ui: React.ReactElement) {
  return render(
    <QueryClientProvider client={queryClient}>
      {ui}
    </QueryClientProvider>
  );
}
```

## Coverage Requirements

**Target: 80% overall coverage**

**Critical paths (require 90%+):**
- `backend/app/services/validation/` - Core validation logic
- `backend/app/api/routes/` - API endpoints
- `frontend/src/hooks/` - React hooks

**Excluded from coverage:**
- Configuration files
- Type definitions
- Generated code

## Test Examples

**Backend - Check test:**
```python
# tests/test_validation/test_basic_checks.py
import pytest
from rdkit import Chem
from app.services.validation.checks.basic import ValenceCheck
from app.schemas.validation import Severity

class TestValenceCheck:
    def setup_method(self):
        self.check = ValenceCheck()

    def test_valid_molecule_passes(self, valid_mol):
        result = self.check.run(valid_mol)

        assert result.passed is True
        assert result.severity == Severity.PASS
        assert result.affected_atoms == []

    def test_invalid_valence_fails(self, invalid_valence_mol):
        result = self.check.run(invalid_valence_mol)

        assert result.passed is False
        assert result.severity == Severity.CRITICAL
        assert len(result.affected_atoms) > 0
```

**Backend - API test:**
```python
# tests/test_api/test_validation_routes.py
import pytest
from httpx import AsyncClient
from app.main import app

@pytest.mark.asyncio
async def test_validate_valid_smiles():
    async with AsyncClient(app=app, base_url="http://test") as client:
        response = await client.post(
            "/api/v1/validate",
            json={"molecule": "CCO", "format": "smiles"}
        )

    assert response.status_code == 200
    data = response.json()
    assert data["overall_score"] == 100
    assert len(data["issues"]) == 0

@pytest.mark.asyncio
async def test_validate_invalid_smiles():
    async with AsyncClient(app=app, base_url="http://test") as client:
        response = await client.post(
            "/api/v1/validate",
            json={"molecule": "invalid_smiles", "format": "smiles"}
        )

    assert response.status_code == 400
```

**Frontend - Component test:**
```typescript
// components/validation/ValidationResults.test.tsx
import { render, screen } from '@testing-library/react';
import { ValidationResults } from './ValidationResults';

const mockResult = {
  overall_score: 85,
  ml_readiness_score: 90,
  molecule_info: {
    canonical_smiles: 'CCO',
    inchikey: 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
  },
  issues: [
    {
      check_name: 'stereo',
      passed: false,
      severity: 'warning',
      message: 'Undefined stereocenter',
    },
  ],
  execution_time_ms: 45,
};

describe('ValidationResults', () => {
  it('displays overall score', () => {
    render(<ValidationResults result={mockResult} />);

    expect(screen.getByText('85')).toBeInTheDocument();
  });

  it('displays issues when present', () => {
    render(<ValidationResults result={mockResult} />);

    expect(screen.getByText('Undefined stereocenter')).toBeInTheDocument();
  });
});
```

---

*Testing analysis: 2026-01-20*

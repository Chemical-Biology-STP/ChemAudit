# Coding Conventions

**Analysis Date:** 2026-01-20

## Code Style

**Python (Backend):**
- Formatter: Black (default settings)
- Import sorter: isort
- Type checker: mypy (strict mode)
- Linter: Configured via pre-commit hooks
- Line length: 88 (Black default)

**TypeScript (Frontend):**
- Formatter: Prettier
- Linter: ESLint with TypeScript rules
- Strict mode enabled in `tsconfig.json`
- Line length: 100

## Naming Patterns

**Python:**
```python
# Classes - PascalCase
class ValidationEngine:
    pass

class UndefinedStereoCentersCheck(BaseCheck):
    pass

# Functions/methods - snake_case
def validate_molecule(mol: Chem.Mol) -> ValidationResult:
    pass

# Variables - snake_case
validation_result = engine.validate(mol)

# Constants - UPPER_SNAKE_CASE
MAX_BATCH_SIZE = 100000
DEFAULT_CHECKS = ["parsability", "valence"]

# Private - leading underscore
def _calculate_score(self, results):
    pass
```

**TypeScript:**
```typescript
// Components - PascalCase
function ValidationResults({ result }: Props) {}

// Hooks - camelCase with use prefix
function useValidation() {}

// Variables/functions - camelCase
const validationResult = await validateMolecule(smiles);

// Types/Interfaces - PascalCase
interface ValidationResponse {}
type Severity = 'critical' | 'error' | 'warning';

// Constants - UPPER_SNAKE_CASE
const API_BASE_URL = '/api/v1';
```

## Code Patterns

**Validation Check Pattern:**
```python
class ExampleCheck(BaseCheck):
    name = "example_check"
    description = "Description of what this check does"
    category = "basic"  # basic, stereo, alerts, representation

    def run(self, mol: Chem.Mol) -> CheckResult:
        # Perform check logic
        passed = True  # or False based on check

        return CheckResult(
            check_name=self.name,
            passed=passed,
            severity=Severity.PASS if passed else Severity.WARNING,
            message="Result message",
            affected_atoms=[]  # List of atom indices if applicable
        )
```

**API Route Pattern:**
```python
from fastapi import APIRouter, HTTPException
from app.schemas.validation import ValidationRequest, ValidationResponse

router = APIRouter()

@router.post("/validate", response_model=ValidationResponse)
async def validate_molecule(request: ValidationRequest):
    try:
        mol = parse_molecule(request.molecule, request.format)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))

    result = validation_engine.validate(mol, checks=request.checks)
    return result
```

**React Hook Pattern:**
```typescript
export function useValidation() {
  const mutation = useMutation({
    mutationFn: (request: ValidationRequest) =>
      validationApi.validate(request),
  });

  return {
    validate: mutation.mutate,
    validateAsync: mutation.mutateAsync,
    isLoading: mutation.isPending,
    error: mutation.error,
    result: mutation.data,
  };
}
```

**React Component Pattern:**
```typescript
interface Props {
  result: ValidationResponse;
}

export function ValidationResults({ result }: Props) {
  // Component logic
  return (
    <div className="space-y-4">
      {/* JSX */}
    </div>
  );
}
```

## Error Handling

**Python:**
```python
# Custom exceptions
class ValidationError(Exception):
    def __init__(self, message: str, details: dict = None):
        self.message = message
        self.details = details or {}

# In route handlers
try:
    result = validate(mol)
except ValidationError as e:
    raise HTTPException(status_code=400, detail=e.message)
except Exception as e:
    logger.exception("Unexpected error")
    raise HTTPException(status_code=500, detail="Internal error")
```

**TypeScript:**
```typescript
try {
  const result = await validationApi.validate(request);
  return result;
} catch (error) {
  if (axios.isAxiosError(error)) {
    throw new Error(error.response?.data?.detail || 'Validation failed');
  }
  throw error;
}
```

## Import Order

**Python (isort configured):**
```python
# 1. Standard library
import json
from pathlib import Path
from typing import List, Optional

# 2. Third-party
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from rdkit import Chem

# 3. Local
from app.core.config import settings
from app.schemas.validation import ValidationRequest
from app.services.validation.engine import validation_engine
```

**TypeScript:**
```typescript
// 1. React/framework
import React, { useState, useEffect } from 'react';

// 2. Third-party
import { useQuery } from '@tanstack/react-query';
import axios from 'axios';

// 3. Local components
import { Button } from '@/components/ui/button';
import { ValidationResults } from './ValidationResults';

// 4. Types
import type { ValidationResponse } from '@/types/validation';

// 5. Styles
import './styles.css';
```

## Documentation

**Python docstrings (Google style):**
```python
def validate_molecule(mol: Chem.Mol, checks: List[str] = None) -> ValidationResult:
    """
    Run validation checks on a molecule.

    Args:
        mol: RDKit molecule object to validate.
        checks: List of check names to run. If None, runs all checks.

    Returns:
        ValidationResult with all check outcomes and scores.

    Raises:
        ValueError: If molecule is None or invalid.
    """
```

**TypeScript (JSDoc):**
```typescript
/**
 * Hook for molecule validation.
 * @returns Object with validate function and state
 */
export function useValidation() {
  // ...
}
```

---

*Conventions analysis: 2026-01-20*

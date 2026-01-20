# Technical Concerns

**Analysis Date:** 2026-01-20

## Project Status

**Status:** Greenfield - Pre-Implementation

This is a new project with comprehensive documentation but no implementation code yet. The concerns listed below are anticipated risks and architectural decisions to address during development.

## Technical Debt

**Current:** None (no code exists yet)

**Anticipated debt from Phase 1 rush:**
- Hardcoded configuration values
- Incomplete error handling
- Missing edge case tests
- Deferred performance optimization

## Known Issues

**None currently** - project is in planning phase.

## Security Considerations

**Input Validation:**
- SMILES strings can contain special characters - sanitize before parsing
- File uploads (SDF, CSV) need size limits and format validation
- Prevent code injection via molecule strings (documented in `docs/ARCHITECTURE.md`)

**Authentication (Phase 3):**
- API key authentication planned
- Rate limiting required (slowapi)
- No user accounts in Phase 1-2 (anonymous access)

**Data Protection:**
- No sensitive data storage in Phase 1
- Uploaded files need secure temp storage with TTL
- Consider GDPR compliance for production

## Performance Risks

**RDKit.js Browser Performance:**
- WebAssembly module is large (~10MB)
- Risk: Slow initial load, memory issues with large molecules
- Mitigation: Server-side rendering fallback, lazy loading

**Batch Processing Scale:**
- Target: 100,000 molecules per batch
- Risk: Memory exhaustion, slow processing
- Mitigation: Chunked processing, Celery workers, progress checkpoints

**Alert Pattern Matching:**
- 1000+ SMARTS patterns to match
- Risk: O(n*m) complexity bottleneck
- Mitigation: RDKit FilterCatalog optimization, caching

**Database Query Performance:**
- Batch results can be large (100K rows)
- Risk: Slow pagination, memory issues
- Mitigation: Proper indexing, streaming responses

## Fragile Areas

**RDKit Version Compatibility:**
- Both Python RDKit and RDKit.js must produce consistent results
- Version mismatches could cause frontend/backend discrepancies
- Pin versions carefully in dependencies

**Stereochemistry Detection:**
- Different RDKit versions handle stereo differently
- Edge cases: undefined stereocenters, pseudo-asymmetric centers
- Requires comprehensive test suite

**Tautomer Canonicalization:**
- Tautomer handling affects InChI/InChIKey generation
- Different canonicalization = different identifiers
- Document and test expected behavior

**External Alert Patterns:**
- PAINS, BRENK, etc. patterns from external sources
- Updates to patterns could change validation results
- Version and track pattern sets

## Architecture Decisions to Validate

**Celery vs. Background Tasks:**
- Decision: Celery for batch processing
- Alternative: FastAPI BackgroundTasks for simpler cases
- Validate: Is Celery complexity justified for expected scale?

**PostgreSQL for Results:**
- Decision: Store all batch results in PostgreSQL
- Alternative: File-based storage for very large batches
- Validate: Query performance at 100K row scale

**Redis for Everything:**
- Decision: Redis for cache, queue, and sessions
- Alternative: Separate services for different concerns
- Validate: Single Redis instance sufficient for load?

## Dependency Risks

**Critical Dependencies:**
- RDKit - Core chemistry, no alternative
- chembl-structure-pipeline - Standardization, could implement in-house if needed
- MolVS - Additional validation, optional

**Version Pinning Required:**
- RDKit (Python): ^2023.9.0
- @rdkit/rdkit (JS): ^2023.9.1
- Keep in sync to avoid discrepancies

**Supply Chain:**
- All dependencies are well-established open source
- RDKit maintained by active community
- No known security advisories for current versions

## Monitoring Gaps (To Address in Phase 4)

**Planned:**
- Prometheus metrics endpoint
- Grafana dashboards
- Sentry error tracking
- Structured logging with correlation IDs

**Key Metrics to Track:**
- Validation latency (P50, P95, P99)
- Batch job throughput
- Error rates by check type
- Queue depth and worker utilization

## Migration Considerations

**Database Schema:**
- Use Alembic for all migrations
- Plan for result schema evolution
- Consider backward compatibility for API changes

**API Versioning:**
- Start with `/api/v1/`
- Plan deprecation strategy before v2

---

*Concerns analysis: 2026-01-20*

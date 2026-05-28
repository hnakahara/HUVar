## Summary

<!-- One or two sentences: what does this PR change? -->

## Type of change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that changes existing TSV output / classification behaviour)
- [ ] Documentation / chore

## Classification impact

If this PR changes evaluator logic, the Bayesian point map, or the 2015
combination rules, briefly describe the expected effect on existing
variants. Example: *"Adds inheritance-aware BS2; ~1% of AR-gene variants
previously called LB may move to Benign."*

If there is **no** classification impact, write "None".

## Test plan

- [ ] `ruff check src tests` passes
- [ ] `pytest tests/unit -v` passes
- [ ] Added unit tests for new logic (or updated existing tests)
- [ ] (If touching annotation logic) ran `pytest tests/integration -v` against
      a populated `data/` directory
- [ ] Updated `CHANGELOG.md` under `## [Unreleased]`

## Linked issues

Closes #

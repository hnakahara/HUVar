# Contributing to HUHVar

Thank you for your interest in improving HUHVar (`acmg-classifier`). This
document describes the development workflow and how to submit changes.

## Reporting issues

Use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.md) for
defects and [feature request template](.github/ISSUE_TEMPLATE/feature_request.md)
for proposed enhancements. **Do not paste real patient data into issues** —
synthesise an equivalent variant from public ClinVar entries instead.

## Development setup

```bash
git clone https://github.com/hnakahara/HUHVar.git
cd HUHVar

mamba create -n acmg-dev -c bioconda -c conda-forge \
    python=3.12 samtools tabix bcftools htslib ensembl-vep=111
mamba activate acmg-dev

pip install -e ".[dev]"
```

## Branching

- `main` — always releasable
- `feat/<short-name>` — new functionality
- `fix/<short-name>` — bug fixes
- `docs/<short-name>` — documentation only

## Pre-commit checklist

Before opening a PR, please make sure:

```bash
ruff check src tests
pytest tests/unit -v
```

both pass cleanly. Integration tests (`tests/integration/`) require a
populated `data/` directory; running them locally is encouraged when your
change touches annotation logic.

## Commit messages

Conventional commits format:

```
<type>(<scope>): <subject>

<optional body explaining WHY, not WHAT>
```

`type` is one of `feat`, `fix`, `refactor`, `docs`, `test`, `chore`, `perf`,
`ci`. Keep the subject under 70 characters.

## Pull requests

- One logical change per PR.
- Reference any related issue (`Closes #NN`).
- Describe the **clinical / interpretation impact** in the PR body if you
  change classifier behaviour. A reviewer needs to understand whether the
  change could move existing variants between categories (e.g. VUS → LP).
- Include a test that fails before your change and passes after, whenever
  practical.

## Adding a new ACMG criterion or strength rule

1. Implement the evaluator under `src/acmg_classifier/criteria/{pathogenic,benign}/`.
2. Register it in the appropriate `CriteriaRegistry` (look at how `PM2`/`BA1`
   are registered).
3. Add a unit test under `tests/unit/test_criteria_*.py` covering at minimum:
   - the positive case (criterion fires at expected strength)
   - the negative case (criterion does NOT fire)
   - one edge case (boundary condition, missing annotation, etc.)
4. Update `README.md` if the change affects user-visible behaviour.

## License of contributions

By submitting a contribution you agree that it will be released under the
project's [Apache License 2.0](LICENSE).

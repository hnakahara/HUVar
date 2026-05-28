---
name: Bug report
about: Report a defect in classification, annotation, or the CLI
title: "[BUG] "
labels: bug
assignees: ''
---

## Summary

A clear, one-sentence description of the bug.

## Environment

- HUHVar (`acmg-classifier`) version (`acmg-classify --version`):
- Python version (`python --version`):
- OS / distribution:
- Assembly (GRCh37 / GRCh38):
- Data source versions (`acmg-classify status`):

## Steps to reproduce

```bash
# Exact command(s) you ran
acmg-classify classify ...
```

## Expected behaviour

What you expected to happen.

## Actual behaviour

What actually happened. Include the relevant **anonymised** TSV row, or the
exact stderr trace.

> Do not paste real patient identifiers. Reproduce the issue with a
> public ClinVar variant or use a synthetic equivalent.

## Additional context

Any other context: gnomAD AF, ClinVar entry, transcript ID, etc.

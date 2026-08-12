# Validation run log

- Run 1 stopped before output writing during the FFPE per-gene audit. The
  legacy helper required a `direction_stability` metadata column that had been
  omitted from the compact validation-panel export. No panel gene, direction,
  weight, size rule or validation result was changed.
- Corrective action: propagate the already-existing discovery
  `direction_stability` field into the helper input, then rerun the unchanged
  validation workflow from the beginning.

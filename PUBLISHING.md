# Publishing the corrected v28.0 public working draft

## GitHub

1. Replace the changed metadata and manuscript files and add
   `RELEASE_NOTES_v28.0.1.md`. Keep the two PNG figures in the repository root.
2. Commit with:

   ```text
   Fix v28.0 citation metadata and figure links
   ```

3. Create release/tag `v28.0-public-draft.1` targeting `main`, with title:

   ```text
   v28.0 public working draft — corrected package
   ```

4. Use `RELEASE_NOTES_v28.0.1.md` as the description, mark the release as a
   pre-release, and attach the corrected full-package ZIP and the manuscript
   PDF.

## Zenodo

GitHub–Zenodo integration is enabled. Publishing the corrective GitHub release
should create another version automatically; do not create a separate manual
record. The record-family DOI is
[10.5281/zenodo.23292159](https://doi.org/10.5281/zenodo.23292159).
The earlier versions, including v10.6 at DOI
[10.5281/zenodo.23292160](https://doi.org/10.5281/zenodo.23292160),
remain available.

Check that the new Zenodo version has version `v28.0-public-draft.1`, preserves
the working-draft description, and includes the source archive containing the
two PNG figures and corrected metadata.

## Final external checks

- The public page calls v28.0 a working draft, not a final paper.
- It does not claim that v10.6 or v13.3 reproduces every v28.0 number.
- The old versions remain accessible.
- No local credentials, tokens, cookies, `.git` directory, or internal audit
  correspondence are uploaded.

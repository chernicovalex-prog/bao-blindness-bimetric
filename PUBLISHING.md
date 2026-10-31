# Publishing the v28.0 public working draft

## GitHub

1. Review `RELEASE_NOTES_v28.0.md`, `README.md`, and `REPRODUCIBILITY.md`.
2. Add the new manuscript, PDF, documentary archive, figures, and metadata to
   the existing repository. Do not delete v21.5, v21.6, or archive v10.6.
3. Commit with a message such as:

   ```text
   Add v28.0 public working draft and version-specific reproducibility notes
   ```

4. Create a new release/tag `v28.0-public-draft`. Suggested title:

   ```text
   v28.0 public working draft — audit-qualified per-probe study
   ```

5. Use `RELEASE_NOTES_v28.0.md` as the release description. Attach the checked
   release ZIP if desired.

## Zenodo

Create a **new version** of the existing record
[10.5281/zenodo.23292160](https://doi.org/10.5281/zenodo.23292160).
Do not replace or erase the published v10.6 version. Copy the metadata from
`ZENODO_METADATA.md`, upload the checked release files, and publish the new
version only after reviewing the file list and description.

If GitHub–Zenodo integration is enabled, a new GitHub release may create the
Zenodo version automatically. Check that the resulting description preserves
the working-draft status and the version-specific reproducibility limits.

## Final external checks

- The public page must call v28.0 a working draft, not a final paper.
- It must not claim that v10.6 or v13.3 reproduces every v28.0 number.
- The old versions must remain accessible.
- No local credentials, tokens, cookies, `.git` directory, or internal audit
  correspondence may be uploaded.

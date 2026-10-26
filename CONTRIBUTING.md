# Branches and releases

| Branch | What it holds | Who changes it |
|---|---|---|
| `main` | Production: what the website serves (GitHub Pages builds from `main`, folder `docs/`) and the current release in `releases/` | Only by merging `develop` once a change is built and tested |
| `develop` | The working branch: every change lands here first | Merges of topic branches |
| `fix/<topic>`, `feature/<topic>` | One change each, branched from `develop` (for example `fix/review-v1.3.1`) | Its author; merged into `develop`, then deleted |
| `release/vX.Y.Z` | The frozen state of each published version: the commit that published `releases/vX.Y.Z/` and its website | Never, once created |

Each published version also has a tag `vX.Y.Z` and a GitHub release (which Zenodo archives). Create the tag from the
GitHub Releases page, targeting the version's `release/vX.Y.Z` branch.

The archive branches point at:

| Version | Commit |
|---|---|
| `release/v1.0.0` | the `v1.0.0` tag (the published run) |
| `release/v1.1.0` | `1b62772`, which published `releases/v1.1.0/` |
| `release/v1.2.0` | `e3ffa38`, which published `releases/v1.2.0/` |
| `release/v1.3.0` | `8371a5a`, v1.3.0 as reviewed and served |
| `release/v1.3.1` | the commit that published `releases/v1.3.1/` (the corrections after the review of v1.3.0) |
| `release/v1.4.0` | the commit that published `releases/v1.4.0/` (the Lu 2024 sites and five annotation layers) |

Older names are kept for their links and may be deleted: `release-v1.1.0` (merged into `main`), `v1.0.0_archive` (the
repository before the reorganisation; the README links it), `inputs-april-2026-update` and the `claude/*` branches.

## Making a change

1. Branch from `develop`: `git switch -c fix/<topic> origin/develop`.
2. Change the code first and commit it, then rebuild what it affects from that clean commit, so the provenance in
   `release.json` and the website data names a commit without local changes.
3. A change to release data is a new version (`v1.3.1` after `v1.3.0`): it gets its own `releases/vX.Y.Z/` with a
   `CHANGES.md`. A published version's folder is never rewritten.
4. Check: `md5sum -c MD5SUMS` inside the release folder, and the website served locally
   (`python -m http.server` inside `docs/`).
5. Merge into `develop`; when the version is ready, merge `develop` into `main`, create `release/vX.Y.Z` from that
   commit, and publish the GitHub release from it.

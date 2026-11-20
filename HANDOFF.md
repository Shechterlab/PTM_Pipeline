# Handoff between PTM_Pipeline and ptmcontext

Two repositories and three roles make a release. This file says who does what, and how commits, the "ok to tag" and the
DOI pass between them. ptmcontext keeps a mirror of it.

## Roles

| Role | Owns | Does not |
|---|---|---|
| **PTM_Pipeline maintainer** (the data side) | `inputs/`, the staging and release scripts in `scripts/`, `releases/`, the version files (`.zenodo.json`, `CITATION.cff`, the README's version row) and their tests | change ptmcontext; merge, tag or publish |
| **ptmcontext maintainer** (the website side) | the ptmcontext repository: the explorers, their build, their data pins and CI | change PTM_Pipeline; merge, tag or publish |
| **Owner** | decisions on scope and sources; every merge into `develop` and `main` in both repositories; the `release/vX.Y.Z` branch, the `vX.Y.Z` tag and the GitHub release | |

Each maintainer works on its own branch and opens pull requests. Only the owner merges.

## A PTM_Pipeline release

1. **Change.** The data side builds the change on a topic branch: staging, `releases/vX.Y.Z/`, data dictionary, tests. It
   opens a pull request into `develop` and sends the website side:
   - the head commit;
   - the new files and their columns;
   - the counts.
2. **Merge into develop.** The owner merges. The data side opens the `develop` → `main` pull request; the owner merges it.
3. **"ok to tag".** The website side builds its explorers from that `main` commit and checks the release files against
   them. It then sends either:
   - "ok to tag", naming the commit; or
   - a concrete change request. That change goes on `develop` and through steps 1–2 again before any tag.
4. **Tag and publish.** The owner then:
   - creates `release/vX.Y.Z` from the `main` commit;
   - tags `vX.Y.Z` on `main` at that commit;
   - publishes the GitHub release.

   Tags stay `vX.Y.Z` on `main`: ptmcontext's CI fails if a site's data pin is not a tagged PTM_Pipeline release.
5. **DOI.** Zenodo archives the tagged commit and mints a version DOI under the all-versions DOI
   (10.5281/zenodo.22696798). The data side then:
   - sends the website side the tag, the commit hash and the version DOI;
   - adds the DOI to the README and `CITATION.cff` in a small pull request into `develop`.
6. **Pin.** The website side pins both the tag name and the commit hash, and cites the version DOI.

## Rules both sides keep

- **Version numbers:**
  - a minor release (v1.7.0) adds sources, marks or layers;
  - a patch release (v1.6.1) corrects or rebuilds without new content.
- **Published releases:** a published release folder is never rewritten. A correction goes into the next release, with a
  note in its `CHANGES.md`.
- **Holds:** nothing is tagged while a change request from the website side is open. An "ok to tag" can be withdrawn
  until the tag exists.
- **Communication:** the two sides talk through their sessions' messages, not files committed to either repository.
  Decisions relayed from the owner are named as such.
- **Data standard:** new data follow the shared standard,
  [`docs/DATA_STANDARD.md`](https://github.com/Shechterlab/ptmcontext/blob/develop/docs/DATA_STANDARD.md) in ptmcontext.
  `ADDING_A_SOURCE.md` is the guide to the steps here.

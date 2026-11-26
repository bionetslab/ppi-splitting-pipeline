## What is this PR about?

<!-- Keep a Changelog (https://keepachangelog.com): one bullet per user-facing change, link issues/PRs.
     Delete empty sections. A release PR collects the entries of the PRs it contains. -->

### Added

### Changed

### Deprecated

### Removed

### Fixed

### Security

## Checklist

- [ ] Pipeline tests pass (or the PR only touches files the tests skip)
- [ ] README updated if parameters, samplesheet columns or outputs changed
- [ ] Dependency changes go into `pixi.toml`, not `environment.yml`; `pixi.lock` is committed alongside
- [ ] `environment.yml` re-exported after any `pixi.toml` change:
      `pixi workspace export conda-environment -p linux-64 -n ppi-splitting-pipeline environment.yml`

## Pre-release (PRs into `main` for a release only)

- [ ] `manifest.version` in `nextflow.config` bumped above the last release (checked by CI)
- [ ] Dependencies refreshed if wanted: `pixi update`, then re-export `environment.yml` (see above)
- [ ] Both images build from the lock: `docker build -f docker/Dockerfile .` and `docker build -f docker/Dockerfile.split_only .`
- [ ] PR title names the version

After merging:

- [ ] Publish a GitHub release whose tag is exactly `manifest.version`
- [ ] `publish-docker.yml` succeeded and both ghcr packages are public

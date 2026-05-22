# Changelog

All notable changes to SMART are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and SMART follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] — 2026-05-22 — First public release

### Added

- First public release on github.com/iterorganization, distributed under
  LGPL-3.0-or-later.
- IMAS actor (`libsmart.a` + iWrap Python wrapper `run_smart.py`) consuming
  `equilibrium`, `core_profiles` and optionally `pellets` IDSs and returning
  an updated `core_profiles` IDS.
- Standalone Fortran executable (`smart`) for IDS-in / IDS-out runs without
  Python.
- Kuteev/Parks ablation model with Strauss-style post-ablation mass
  relocation along the magnetic field, following Polevoi & Shimada
  (PPCF **43**, 2001, 1525).
- Configurable ECRH source term (`ECH2a` parameter block).
- XSD-validated XML input files for both standalone and actor entry points.
- ITER SDCC build scripts under `ci-sdcc/`.
- `CITATION.cff`, `CONTRIBUTING.md`, `COPYING`, `COPYING.LESSER`,
  per-source-file SPDX headers.

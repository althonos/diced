# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/diced/compare/v0.1.2...HEAD


## [v0.1.2] - 2024-11-04
[v0.1.2]: https://github.com/althonos/diced/compare/v0.1.1...v0.1.2

### Changed
- Bump `pyo3` dependency to `v0.22.5`.
- Use `maturin` to build the Python package.
- Use PyData theme for the Sphinx documentation of the Python package.


## [v0.1.1] - 2024-06-19
[v0.1.1]: https://github.com/althonos/diced/compare/v0.1.0...v0.1.1

### Changed
- Use raw sequence bytes to avoid panics on Unicode characters slicing in `Scanner`.

### Fixed
- Incorrect metadata in Python package and documentation.
- Missing test data in Python wheels.
- Indexing and underflow errors found by `afl` fuzzer.


## [v0.1.0] - 2024-06-11
[Unreleased]: https://github.com/althonos/diced/compare/11ad0d3...v0.1.0

Initial release.

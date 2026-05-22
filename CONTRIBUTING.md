# Contributing to SMART

Thanks for your interest in SMART. This is a small physics code maintained by
the ITER Organization Plasma Modelling & Analysis section; contributions of
all kinds — bug reports, documentation fixes, performance improvements, new
features — are welcome.

## Reporting bugs and asking questions

Please open an issue on
[github.com/iterorganization/smart](https://github.com/iterorganization/smart/issues)
with:

- a clear description of what you observed and what you expected,
- the SMART version (commit hash or release tag),
- the compiler and IMAS Access Layer versions you are using,
- a minimal scenario (IDS pulse/run, or the relevant `input/*.xml`) that
  reproduces the problem, where possible.

## Submitting changes

1. Fork the repository and create a topic branch from `main`.
2. Make your change. Keep commits focused; prefer one logical change per
   commit.
3. Run `make validate` and rebuild (`make clean && make`) to make sure both
   the Fortran library and the standalone executable still compile cleanly
   with your toolchain.
4. Add or update tests where applicable.
5. Open a pull request against `main`. Describe **what** changed and
   **why**; link to any related issue.

## Developer Certificate of Origin

By contributing to SMART you certify that the contribution is your own work
(or that you have the right to submit it under the project licence) and that
the contribution is licensed under LGPL-3.0-or-later — i.e. you agree to the
[Developer Certificate of Origin v1.1](https://developercertificate.org/).
Please sign off each commit with `git commit -s`, which appends a
`Signed-off-by:` line.

## Coding style

- **Fortran**: prefer free-form Fortran 90 for new files; match the
  surrounding style when editing fixed-form `.f` files. `implicit none` is
  required in new code.
- **Python**: PEP 8 where reasonable; avoid pulling in heavy new
  dependencies.

## License

By contributing, you agree that your contributions will be licensed under the
GNU Lesser General Public License, version 3.0 or later (LGPL-3.0-or-later).

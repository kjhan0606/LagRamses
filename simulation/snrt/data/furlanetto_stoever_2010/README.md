# Furlanetto--Stoever 2010 fast-electron tables

These are the 14 electronic tables used for the Furlanetto & Stoever (2010)
fast-electron energy-deposition calculation:

- S. R. Furlanetto & S. J. Stoever, *Secondary ionization and heating by fast
  electrons*, MNRAS 404, 1869 (2010),
  [doi:10.1111/j.1365-2966.2010.16401.x](https://doi.org/10.1111/j.1365-2966.2010.16401.x),
  [arXiv:0910.4410](https://arxiv.org/abs/0910.4410).

The files were taken from the official
[`21cmFAST`](https://github.com/21cmfast/21cmFAST) distribution at commit
`892f98c80cfe985ca6b399ec6b51a3aa95124b11`, under
`src/py21cmfast/_data/x_int_tables`. Its interpolation source identifies
Steven Furlanetto as the author. The data cover 258 electron energies from
10 to 9937.21 eV and 14 H II fractions from `1e-4` to `0.999`.

The numeric/text payload of every upstream file is unchanged. A single final
LF was appended because the upstream files omit a terminal newline. Both
upstream and vendored SHA256 values are recorded in `TABLE_MANIFEST.json` and
checked by the SNRT artifact test.

The calculation assumes primordial H/He abundance, equal H II and He II
fractions, and negligible He III. SNRT uses the actual H II fraction as the
table coordinate, preserves H I, He I, and He II secondary-ionization channels
separately, and records this composition approximation explicitly.

The upstream files are distributed under the 21cmFAST MIT license reproduced
in `LICENSE-21cmFAST.txt`.

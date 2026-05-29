[1]: https://ramses-organisation.readthedocs.io/en/latest
[5]: https://bitbucket.org/vperret/dice
[6]: https://bitbucket.org/ohahn/music
[7]: https://github.com/osyris-project/osyris
[8]: https://github.com/pynbody/pynbody
[9]: https://yt-project.org

# cuRAMSES

**cuRAMSES** is an in-house, performance-optimised fork of the
[RAMSES](https://bitbucket.org/rteyssie/ramses) cosmological AMR code,
developed at the Korea Institute for Advanced Study (KIAS) for the
Horizon Run 5 and HR6 simulations. The cuRAMSES contributions
(recursive *k*-section domain decomposition, Morton-key hash-table
neighbour lookup, FFTW3 direct Poisson solver, spatial hash-binning
for SN/AGN feedback, variable-rank HDF5 restart, and CUDA-accelerated
multigrid/Godunov kernels) are bundled under [`patch/cuRamses/`](./patch/cuRamses/).

### Upstream lineage

cuRAMSES is derived from
[`bitbucket.org/rteyssie/ramses`](https://bitbucket.org/rteyssie/ramses)
master branch at commit
[`d92cd7d9`](https://bitbucket.org/rteyssie/ramses/commits/d92cd7d9)
(2016-11-08, RAMSES version 3.10), as identified by `git`-blob-hash
intersection of all unmodified base-directory source files. The
sub-grid physics (sink particles, AGN feedback, supernova feedback,
metal yields) trace back to the Horizon-AGN version of RAMSES
maintained by Yohan Dubois and collaborators, which is itself
not publicly distributed.

### Companion paper

J. Kim, *cuRAMSES: Scalable AMR Optimizations for Large-Scale
Cosmological Simulations*, MNRAS (2026, in revision).
Cite this paper when using cuRAMSES.

### Acknowledgements

The cuRAMSES contributions in [`patch/cuRamses/`](./patch/cuRamses/)
rest on contributions from many people:

- **Romain Teyssier** — public release and continued development of RAMSES.
- **Yohan Dubois** — Horizon-AGN version of RAMSES, on which the present
  in-house lineage is based. We thank him for sharing his code.
- **Taysun Kimm** — chemical enrichment model (alpha-element and stellar-wind yields).
- **Christophe Pichon, Sébastien Peirani, and Julien Devriendt** — design
  discussions on the Horizon-AGN code and simulation.
- The upstream RAMSES developer community for continued maintenance and
  improvements (see [`CONTRIBUTORS.md`](./CONTRIBUTORS.md) for the
  full list).

cuRAMSES is supported by the Korea Institute for Advanced Study and the
KIAS Center for Advanced Computation. J.~Kim is supported by KIAS
Individual Grants (KG039603).

### Build

A reference build is provided as [`bin/Makefile`](./bin/Makefile) /
[`patch/cuRamses/Makefile`](./patch/cuRamses/Makefile). The default
configuration uses the Intel `ifx` compiler with HDF5 parallel I/O and
optional CUDA acceleration (`make HDF5=1 USE_CUDA=1`). 128-bit Morton
keys for deep-AMR runs are enabled with `MORTON128=1`.

---

![GitHub logo dark-mode-only](./doc/img/full_project_logo_dark.svg#gh-dark-mode-only)
![GitHub logo light-mode-only](./doc/img/full_project_logo.svg#gh-light-mode-only)

## This is the [ramses](https://github.com/ramses-organisation/ramses/) GitHub repository.

Ramses is an open source code to model astrophysical systems, featuring self-gravitating, magnetised, compressible, radiative fluid flows. It is based  on the Adaptive Mesh Refinement (AMR)  technique on a  fully-threaded graded octree.
[ramses](https://github.com/ramses-organisation/ramses/) is written in  Fortran 90 and is making intensive use of the Message Passing Interface (MPI) library.

### ⬇️ Get the code

Download the code by cloning the git repository using
```
$ git clone https://github.com/ramses-organisation/ramses
```

### ℹ️ Get Support

You can go to the user's guide using [Read The Docs][1].
Please register also to the [mailing list](http://groups.google.com/group/ramses_users).
You can get support on Ramses on the [Github Discussion page](https://github.com/ramses-organisation/ramses/discussions)

### 🛠️ Tools and ressources

To generate idealised initial conditions of galaxies, check out [DICE][5].
To generate cosmological initial conditions, check out [MUSIC][6].

To visualize RAMSES data, we encourage you to use [YT][9], [OSYRIS][7] or [PYNBODY][8].

You'll find a lot of useful ressources, links and news about the Ramses community on https://ramses.cnrs.fr/, a website edited by the French "SNO" Ramses.

# Third-party source packages

This directory is a collection of independently licensed sources. Including a
dependency does not transfer its authorship, and no repository-level notice
relicenses third-party code. Preserve each source package's `DESCRIPTION`,
authorship and license notices when redistributing it.

## Supplied comparator dependencies

| Source directory | Original supplied archive | Declared license | Authors recorded in DESCRIPTION |
| --- | --- | --- | --- |
| `vendor/HMP` | `HMP_2.0.tar.gz` | Apache License (== 2.0) | Patricio S. La Rosa, Elena Deych, Sharina Carter, Berkley Shands, Dake Yang and William D. Shannon |
| `vendor/miLineage` | `miLineage_4.0.tar.gz` | GPL (>= 2) | Zheng-Zheng Tang |

The source directories were extracted from the locally supplied source
archives. R source, help files, namespace declarations and package-bundled
datasets are retained. No fitting functions or comparator algorithms were
modified. The original archives contain no separate `LICENSE` or `NOTICE`
files. For convenience, the standard license texts distributed with R are
included under `vendor/LICENSES`; these are not additional author notices.

The original supplied source archives were identified by these SHA-256 values:

```text
HMP_2.0.tar.gz       4b92388f43c37a6a50903daf6a1ef6c3f7b79ac1cb024972c9963a72abeefe96
miLineage_4.0.tar.gz 88e2ff79e0db9258350a68608f2b6201e58f36f6dbdc082d9c503d994cff1d47
```

HMP's original metadata identifies CRAN as its repository:
<https://cran.r-project.org/package=HMP>. miLineage's supplied metadata
identifies version 4.0, dated 2024-05-09, and its original author/maintainer.
The copy does not claim a newer public release or substitute another package
with a similar name.

## MetaQCAT

MetaQCAT 0.2.1 retains its GPL (>= 3) declaration and `LICENSE.md`.
It is separate from the supplied comparator packages above.

Other dependencies are installed separately from their normal repositories;
their source code is not included merely because an analysis script calls
them. Their own licenses continue to apply.

## Frozen analysis inputs

The frozen source tables under `real/CRC/data/`, `real/IBD/data/` and
`simulation/data/` are distributed as research inputs, not as newly generated
data or association results. Their source and transformations are documented
in `INPUTS.md`. They are not relicensed by the MetaQCAT package license.

CRC profiles and metadata are attributed to Wirbel et al. (2019),
[Zenodo 3517209, version 1.0](https://doi.org/10.5281/zenodo.3517209) and the
[CRC meta-analysis repository](https://github.com/zellerlab/crc_meta).
The bundled taxonomic annotation is the frozen project mOTU-to-genus mapping.

IBD tables are extracted without extra processing from Ma et al. (2022)'s
[IBD analysis repository](https://github.com/biobakery/ibd_paper), fixed at
commit `ae9d98ad170546f02af3c67797cdba971bc95e7d`.
The upstream repository's MIT notice is retained below:

```text
MIT License

Copyright (c) 2022 Siyuan Ma

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

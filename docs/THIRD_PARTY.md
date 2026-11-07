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

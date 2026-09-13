# Vignette calculations

Three vignettes retain precomputed model or bootstrap results so package
checks do not repeat the same calculations. Their displayed code, sample
sizes, seeds and replication counts are unchanged. Ordinary tables and
figures are rebuilt from those results.

From the package root, regenerate the results with:

```sh
Rscript tools/vignettes/precompute.R
Rscript tools/vignettes/precompute.R --check
```

The script executes the vignette code itself, then records the package
version, R version, execution time, source hashes and result-file hashes.
Source hashes use LF line endings on every platform. A source build refuses
outdated results; CI checks the manifest before building. To render a
vignette with fresh calculations, set `RASCH_REBUILD_VIGNETTES=true`.

The recorded fits use R serialization format 3 (R 3.5.0 or later). Their
runtime signatures depend on the R build. The reader verifies the complete
record's checksum before refreshing those signatures for the current R
session; it does not change the saved files or any numerical results.
The generator checks fit/result pairing before writing the manifest.

This does not replace the package tests or the simulation studies in
`tools/simval/`. The full test suite remains enabled in CI with
`NOT_CRAN=true`.

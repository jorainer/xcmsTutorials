# Changelog

## xcmsTutorials 1.2

### Changes in 1.2.0

- Update to Bioconductor release 3.22 and add *xcms* reference.

## xcmsTutorials 1.1

### Changes in 1.1.2

- Add citation information to the README and the package.

### Changes in 1.1.1

- Use
  [`spectraSampleIndex()`](https://rdrr.io/pkg/MsExperiment/man/MsExperiment.html)
  instead of the old
  [`fromFile()`](https://lgatto.github.io/MSnbase/reference/Spectrum-class.html)
  function to get sample assignments of individual spectra.
- Add `()` to all functions in the descriptive text to discriminate
  functions from classes or parameter names.

### Changes in 1.1.0

- Release the tutorial to Zenodo:
  <https://doi.org/10.5281/zenodo.11185521>
- Use Bioconductor 3.19 release docker image and packages.
- Add examples for simple quality assessment of BPC or BPS data.
- Introduce parameter `ppm` for `PeakDensityParam` correspondence
  analysis method.
- Add an additional section for extraction of chromatographic or spectra
  data for identified chromatographic peaks and perform isotopologue
  search in these.

## xcmsTutorials 1.0

### Changes in 1.0.4

- Use the Bioconductor 3.18 release docker image.

### Changes in 1.0.3

- Add reference to the `MsQuality` paper.

### Changes in 1.0.2

- Add the `msdata` to *depends* to ensure it’s going to be installed.

### Changes in 1.0.1

- Small fixes with bullet points in the README and
  corrections/reformulations in the workshop Rmd.

### Changes in 1.0.0

- Restructure and clarify basic data access section.

## xcmsTutorials 0.99

### Changes in xcmsTutorials 0.99.3

- Finalize the additional visualization paragraph.

### Changes in xcmsTutorials 0.99.2

- Add additional visualization options. Require `xcms` version \>=
  3.99.4.

### Changes in xcmsTutorials 0.99.1

- Add requirements for `xcms` and `Spectra` packages as well as R
  version.

## xcmsTutorials 0.2

### Changes in xcmsTutorials 0.2.0

- Rewrite parts of the tutorial based on feedback and corrections
  provided from Philippine Louail.

## xcmsTutorials 0.1

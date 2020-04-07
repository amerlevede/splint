# splint

Splint is an algorithm designed to infer copy number variations in yeast genomes with a bias in read depth. In particular, when the read depth changes continuously across each chromosome (e.g. being lowest in the center and highest near the sides).

The algorithm fits the bias using a smoothing spline non-linear regression, adding extra discontinuous step functions at points where changes in read depth have been inferred. These discontinuities are detected by finding peaks in the differences of rolling-frame-average read depths. The size of the jumps comes from the regression in the form of the step functions of the coefficients.

More information in thesis.pdf.

This algorithm was used for CNV inference in [Gallone et al. Domestication and Divergence of Saccharomyces cerevisiae Beer Yeasts. 2016][https://www.ncbi.nlm.nih.gov/pubmed/27610566].

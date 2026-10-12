# tidyseq
Convenience functions for transcriptomic analyses 

## This version

This branch holds tidyseq 0.1, the version used for the analyses in the AJRCCM ARDS phenotypes paper. It is kept unchanged so that code still runs. Install it with:

```r
pak::pak("AartikSarma/tidyseq@v0.1.2")
```

The `main` branch holds an unrelated rewrite with different functions.

Many of these functions are wrappers around Bioconductor functions. Functions are designed to try to follow the unifying principles described here: https://design.tidyverse.org/unifying.html

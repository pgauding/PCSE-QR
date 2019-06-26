# PCSE_sparse

A repository for final project for POLS 904 at University of Kansas,
FS 17, Paul Johnson professor.

This project demonstrates the feasability of using sparse matrices to
improve management and performance in R. Specifically, I modify and
extend the capabilities of the {pcse} package by implementing sparse
capabilities provided by the {Matrix} package. I identify a recurring
issue in processing large _N x T_ matrices in which the package fails
to produced the required Kronecker multiplication and crashes R. Then,
I identify the matrix storage structure used in the original package,
which I hypothesize causes memory overload. I provide background on
the sparse matrxi alternative, and provide examples of
implementation. I find notable improvements in both dimensional size
and speed capabilities.

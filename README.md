# Multivariate-Functional-PCA
R script for implementing the Multivariate Functional PCA (MFPCA) proposal from Berrendero, Justel, and Svarc in the paper "Principal components for multivariate functional data".
The main function is pwMFPCA() an it is assumed that the discretized points from the functional observations are equidistant.

## Arguments: 
- mfData: Multivariate Functional Data Object in format mfData from *roahd package* available from CRAN.
- ncomp: Number of components for the MFPCA.

## Values:
- functions: A roahd object with the new curves represented by ncomp components.
- pi_1: The cumulative variability by components according to criterion $\pi_1$.
- pi_2: The cumulative variability by components according to criterion $\pi_2$.
- eigen: A list containing the eigenvalues and eigenvectors of the Sigma(t) matrix.

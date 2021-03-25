# highEPs
'highEPs' is a repository to study asymptotic exceptional points in large systems of subwavelength resonators.

Please cite the following reference when using this source:

[1] Ammari H, Davies B, Hiltunen E O, Lee H & Yu s. 2021 *High-order exceptional points and enhanced sensing in subwavelength resonator arrays* Studies in Applied Mathematics 146: 440--462 (https://doi.org/10.1111/sapm.12349)

# highEPs/dilute

This code studies the dilute weighted capacitance matrix, using the approximation for the capacitance coefficients (Lemma 1 in [1]) given by:
```math
\begin{equation}
C_{ij} = \begin{cases}
\mathrm{Cap}_B + O(\epsilon^2), &\quad i=j,
\\
-\epsilon \frac{(\mathrm{Cap}_B)^2}{4\pi |z_i-z_j|} + O(\epsilon^2), &\quad i\neq j,
\end{cases}
\end{equation}
```
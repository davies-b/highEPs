# highEPs
'<highEPs>' is a repository to study asymptotic exceptional points in large systems of subwavelength resonators.

Please cite the following reference when using this source:

[1] Ammari H, Davies B, Hiltunen E O, Lee H & Yu s. 2021 *High-order exceptional points and enhanced sensing in subwavelength resonator arrays.* Studies in Applied Mathematics 146: 440--462 (https://doi.org/10.1111/sapm.12349)

## highEPs/dilute

This code studies the dilute weighted capacitance matrix, using the approximation for the capacitance coefficients (Lemma 1 in [1]) given by:

<img src="https://latex.codecogs.com/svg.latex?\large&space;C_{ij}%20=%20\begin{cases}\mathrm{Cap}_B%20+%20O(\epsilon^2),%20&\quad%20i=j,\\-\epsilon%20\frac{(\mathrm{Cap}_B)^2}{4\pi%20|z_i-z_j|}%20+%20O(\epsilon^2),%20&\quad%20i\neq%20j.\end{cases}">

'<highEPs/dilute/Cd1v.m>' generates the dilute capacitance matrix.

'<highEPs/dilute/RUN_EPdilute.m>' searches for parameter values that give an exceptional point of the dilute capacitance matrix. This was used to generate Figures 2 and 5 in [1].

'<highEPs/dilute/RUN_EP_tau.m>' computes how the eigenvalues of the dilute capacitance matrix vary in a region of the exceptional point (which occurs when <img src="https://latex.codecogs.com/svg.latex?\large&space;\tau=1">). This was used to generate Figure 3 in [1].

## highEPs/multipole

'<highEPs/multipole/RUN_EP_tau.m>' computes how the resonant frequencies vary in a region of the asymptotic exceptional point (which occurs when <img src="https://latex.codecogs.com/svg.latex?\large&space;\tau=1">). This uses the multipole expansion method and was used to generate Figure 4 in [1].

'<highEPs/multipole/RUN_EP_modeplot.m>' computes the eigenmodes of the full differential system at the exceptional point, using the multipole expansion method. This was used to generate Figure 6 in [1].

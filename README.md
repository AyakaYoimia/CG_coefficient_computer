# CG_coefficient_computer
A tool to compute Clebsch-Gordan coefficients.

## Remaining Problem
- All Clebsch-Gordan coefficients are $\sqrt{\frac{p}{q}}$, but the current program computes them as `float`.
- ~~If `J != j_1 + j_2`, the recursion runs error because `j`s are `float`~~

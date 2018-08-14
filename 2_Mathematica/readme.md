## Mathematica code that generates Fortran code

- The two mathematica notebooks here ([element formulation](element_formulation.pdf) and [numerical integration](numerical_integration.pdf)) generate some bulk Fortran code that is needed for the user element implementation. The files are commented and sufficiently self-explanatory.
  * `element_formulation.n` exports the element matrices to determine the polynomial coefficients from the DOF and to to evaluate the shape functions, see the file `code56.txt`
  * `numerical_integration.nb` tests the numerical integration scheme and exports the integration point coordiantes and weights, see the file `29GP.txt`
- The content of the `txt` files needs to be copy-pasted to the user element Fortran file.

Continue here: 3. [Fortran UEL](../3_Fortran_UEL)

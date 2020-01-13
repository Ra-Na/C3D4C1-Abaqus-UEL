# Incompatible mode C<sub>1</sub>-continuous tetrahedral 4-node 3D user element for Abaqus

This repository contains all the necessary tools to use a tetrahedral incompatible mode C<sub>1</sub>-continuous user element in Abaqus that has the following features:

1. Fully compatible with C3D4-meshes
2. Large strain formulation
3. Simple user material interface
4. 4 degrees of freedom (DOF, field value and gradient) at each corner node, i.e.~16 total DOF per scalar field
5. Full 5<sup>th</sup> order polynomial shape function
6. Numerical integration exact up to a full 6<sup>th</sup> order polynomial

It supplements the article **_A C1 incompatible mode element formulation for strain gradient elasticity_** [DOI: 10.1007/978-3-030-30406-5_6](https://doi.org/10.1007/978-3-030-30406-5_6)

## Walkthrough

1. [Element formulation](./1_Element_formulation)
2. [Mathematica: Shape function derivation and code generation](./2_Mathematica)
3. [Fortran UEL](./3_Fortran_UEL)
4. [Abaqus example](./4_ABQ_Example)
5. [Visualisation with Paraview](./5_Visualisation)


## Referencing

Please cite as follows:

```
@article{GLUEGE20192,
  author = {Gl{\"u}ge, R.},
  editor = {Altenbach, H. and M{\"u}ller, W.H. and Abali, B.E.},
  title = {A C1 Incompatible Mode Element Formulation for Strain Gradient Elasticity},
  booktitle = {Higher Gradient Materials and Related Generalized Continua},
  year = {2019},
  publisher = {Springer International Publishing},
  address = {Cham},
  pages = {95--120},
  isbn = {978-3-030-30406-5},
  doi = {10.1007/978-3-030-30406-5_6},
  url = {https://doi.org/10.1007/978-3-030-30406-5_6}
}
```

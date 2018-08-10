# Incompatible mode C<sub>1</sub>-continuous tetrahedral 4-node 3D user element for Abaqus

This repository contains all the necessary tools to use a tetrahedral incompatible mode C<sub>1</sub>-continuous user element in Abaqus that has the following features:

1. Fully compatible with C3D4-meshes
2. Large strain formulation
3. Simple user material interface
4. Allows for numeric determination of element stiffness by applying finite differences to material subroutine (no more DDSDDE)
5. 4 degrees of freedom (DOF, field value and gradient) at each corner node, i.e.~16 total DOF
6. Full 5<sup>th</sup> order polynomial shape function
7. Numerical integration exact up to a full 6<sup>th</sup> order polynomial

## Walkthrough

1. [Element formulation](./1_Element_formulation)
2. [Mathematica: Shape function derivation and code generation](./2_Mathematica)
3. [Fortran UEL](./3_Fortran_UEL)
4. [Abaqus example](./4_ABQ_Example)
5. [Visualisation with Paraview](./5_Visualisation)

## Element formulation

- A C<sub>1</sub>-continuous tetrahedral element requires a complete 9<sup>th</sup>-order polynomial, which has in 3D binomial(12 over 9)=220 monomials, see Ženı́šek, A. (1973). “Polynomial approximation on tetrahedrons in the finite element method”. In: Journal of Approximation Theory 7.4, pp. 334–351.
- Taking the field value &Phi; and its gradient &nabla;&Phi; as 4 DOF, this requires 55 nodes. This is too much. 
- Therefore, we consider an incompatible mode element, requiring C<sub>1</sub>-continuity only at specific nodes
- To be compatible with ordinary C3D4 meshes we extrapolate pseudo-DOF from the four corner nodes 

![Element formulation sketch](./element_formulation.gif "Element formulation sketch")

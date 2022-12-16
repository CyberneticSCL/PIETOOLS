# PIETOOLS
PIETOOLS is a free MATLAB toolbox for manipulating Partial Integral (PI) operators, simulating PDE and Partial Integral Equation (PIE) solutions, 
and solving Linear PI Inequalities (LPIs) which are convex optimization problems involving PI variables and PI constraints. 

PIETOOLS can be used to:
- define PI operators acting on functions in 1 and 2 spatial variables (1D and 2D)
- declare 1D and 2D PI operators variables (postive semidefinite or indefinite)
- add operator inequality constraints
- solve LPI optimization problems
- numerically simulate 1D PDE, DDE, and PIE solutions.
	
The interface is inspired by YALMIP and the program structure is based on that used by SOSTOOLS. By default the LPIs are solved using SeDuMi, however, the toolbox also supports use of other SDP solvers such as Mosek, sdpt3 and sdpnal.

For more information on installation and use refer to the documnetation in PIETOOLS/PIETOOLS_documentation/, and demos in PIETOOLS/PIETOOLS_demos/.

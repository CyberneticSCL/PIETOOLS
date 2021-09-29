# PIETOOLS
PIETOOLS is a free MATLAB toolbox for manipulating Partial Integral (PI) operators and solving Linear PI  Inequalities (LPIs) which are convex optimization problems involving PI variables and PI constraints. 

PIETOOLS can be used to:
- define 3-PI or 4-PI operators
- declare 3-PI or 4-PI operators variables (postive semidefinite or indefinite)
- add operator inequality constraints
- solve LPI optimization problems
	
The interface is inspired by YALMIP and the program structure is based on that used by SOSTOOLS. By default the LPIs are solved using SeDuMi, however, the toolbox also supports use of other SDP solvers such as Mosek, sdpt3 and sdpnal.

For more information on installation and use refer to the files in PIETOOLS/PIETOOLS_examples/ folder path.

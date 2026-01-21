function PDE_out = nonlin_prod(PDE_1,PDE_2)
% PDE_out = nonlin_prod(PDE_1,PDE_2) builds a PDE representing the prduct of the 
% terms in the PDE objects PDE_1 and PDE_2.
%
% INPUT
% - PDE_1:      'pde_struct' object corresponding to either a single PDE
%               variable (state, input, or output) of size n or
%               corresponding to one or several terms in a PDE, belonging
%               to n equations.
% - PDE_2:      'pde_struct' object representing a free term to multiply 
%                the equation specified by PDE_1 by. 
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing either the same equation
%               as in PDE_1 but multiplied by ther term in PDE_2, or a
%               new set free term (collected in PDE_out.free)
%               corresponding to the product of PDE_1 and PDE_2.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% CR, 01/19/2026: Initial coding

PDE_out = PDE_1*PDE_2.free;


end
function Pt = ctranspose(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pt = ctranspose(P) returns the conjuagate transpose of 'nopvar' object P.
% Version: 1.0
% 
% INPUT
% P:  nopvar object;
% 
% OUTPUT
% Pt: nopvar object representing adjoint of P;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding CR - 1/16/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isa(P,'nopvar')
    error('Input must be an opvar variable.');
end

Pt = nopvar();
Pt.C = cell(3,1);

for ii=1:numel(P.C)
    Pt.C{ii} = P.C{ii}.';
end

Pt.deg = P.deg;
 
Pt.vars = [P.vars(:,2) P.vars(:,1)];

Pt.dom = P.dom;



end
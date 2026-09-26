function [d,sp] = cx_hinf_dimsp(dim,PIE)
% [D,SP] = CX_HINF_DIMSP(DIM,PIE) turns a 1-D opvar dim pair [out in]
% (rows R^n, L_2[s]) into the dims/spaces structs lpivar_cdopvar takes for
% a rectangular operator, e.g. the synthesis variable Z = lpivar(prog,
% Buop.dim(:,[2,1]),ddZ). Zero-dimension spaces are dropped, since
% lpivar_cdopvar rejects zero counts. Used by cx_Hinf_control and
% cx_Hinf_estimator.
%
% Initial coding MMP, 09/25/2026
L2 = {cell(1,0), {char(PIE.vars(1,1).varname{1})}};
d = struct();   sp = struct();          % fields assigned one by one:
ko = dim(:,1)>0;    ki = dim(:,2)>0;    % struct() with cells builds arrays
d.out = dim(ko,1);  d.in = dim(ki,2);
sp.out = L2(ko);    sp.in = L2(ki);
end

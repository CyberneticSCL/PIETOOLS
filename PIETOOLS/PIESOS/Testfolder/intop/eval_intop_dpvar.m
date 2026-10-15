function Kop_num = eval_intop_dpvar(Kop,d0map)
% KOP_NUM = EVAL_INTOP_DPVAR(KOP,D0MAP) evaluates the coefficient matrices
% in the 'intop' object KOP at fixed numeric decision-variable values,
% returning an 'intop' without decision variables. Used as ground truth for
% testing 'intop' objects with decision-variable kernels.
%
% INPUTS
% - Kop:    m x n 'intop' object, possibly with decision variables Kop.Zd;
% - d0map:  containers.Map from decision-variable name (as in Kop.Zd) to a
%           numeric value; see also eval_dpvar_at.m;
%
% OUTPUTS
% - Kop_num: m x n 'intop' object, equal to Kop except that the decision
%           variables are replaced by the values in d0map, so that
%           Kop_num.Zd is empty;
%
% NOTES
% - Kop.params{k} holds, per matrix row, the constant term and the
%   coefficients of Kop.Zd, so the values are substituted as
%   kron(I_m,[1,d'])*params{k}, following the definition of the class.

% DJ (CC), 10/15/2026: Initial coding;

Kop_num = Kop;
nd = numel(Kop.Zd);
d0 = zeros(1,nd);
for l=1:nd
    d0(l) = d0map(Kop.Zd{l});
end
sel = kron(speye(size(Kop,1)),[1,d0]);
Kop_num.params = cellfun(@(a) sel*a, Kop.params, 'UniformOutput',false);
Kop_num.Zd = cell(0,1);
Kop_num.chkval = true;

end

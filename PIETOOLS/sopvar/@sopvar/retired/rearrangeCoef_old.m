% RETIRED 2026-08-30. Formerly @sopvar/mtimes.m lines 1016-1105.
% Had no callers anywhere in the repo before removal.

function Cout = rearrangeCoef_old(Cmul,C,p,q,Zint)
% rearrangeCoef
%
% Converts
%
%   (I_p \otimes Zmix(s,t)') * Cmul*C
%
% where
%
%   Zmix(s,t) = kron_i (Z_i(s_i) \otimes Z_i(t_i))
%
% into
%
%   (I_p \otimes Z(s)') * Cout * (I_q \otimes Z(t)).
%
% INPUTS
%   Cmul : typically kron(I_p,Csep')
%   C    : size p*NG by q
%   p    : left matrix dimension
%   q    : right matrix dimension
%   Zint : 1-by-n cell basis for both s and t
%
% OUTPUT
%   Cout : size p*NZ by q*NZ
%
% where NZ = prod_i numel(Zint{i}).

Zint = Zint(:).';
ns = cellfun(@numel,Zint);

if isempty(ns)
    Cout = Cmul*C;
    return
end

NZ = prod(ns);

D = full(Cmul*C);

if ~isequal(size(D),[p*NZ^2,q])
    error('rearrangeCoef: input dimensions are inconsistent.');
end

n = numel(ns);

% D has rows ordered as:
%
%   [mixed monomial index, p-index]
%
% where the mixed monomial basis is
%
%   kron_i (Z_i(s_i) \otimes Z_i(t_i)).
%
% For each variable i, the local mixed basis is ordered as
%
%   kron(Z_i(s_i),Z_i(t_i)),
%
% so its local dimensions are [nt_i, ns_i] in MATLAB reshape order.
%
% First reshape into:
%
%   t_n, s_n, ..., t_1, s_1, p, q
%
% because MATLAB / kron ordering reverses dimensions under reshape.

D = reshape(D,[repelem(fliplr(ns),2), p, q]);

% Separate all s-indices from all t-indices.
%
% Current dimensions:
%   t_n, s_n, t_{n-1}, s_{n-1}, ..., t_1, s_1, p, q
%
% Desired dimensions:
%   s_n, ..., s_1, p, t_n, ..., t_1, q
%
% This gives row block (s,p) and column block (t,q).

perm = [2:2:2*n, 2*n+1, 1:2:2*n, 2*n+2];

D = permute(D,perm);

Cout = reshape(D,[NZ*p, NZ*q]);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

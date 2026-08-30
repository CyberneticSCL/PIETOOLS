function [C_gam_alp_beta,ZL,ZR] = int_semisep(G,idxbeta,idxalpha,lims,Csize)
% int_semisep
%
% Evaluates the two-sided semiseparable integral of Sec. 6.1 of the sopvar
% document,
%
%   int_{s''} I_beta(s-s'') I_alpha(s''-s') G(s'') ds''
%       = sum_gamma (Delta_{p(gamma,alpha,beta)} G)(s,s') I_gamma(s-s')
%
% converting
%
%   G(s3a) = (I_g1 \otimes ZG(s3a)')*CG
%
% into semiseparable coefficient matrices
%
%   I_{gamma,beta,alpha}G
%      = (I_g1 \otimes ZL(s3a)') ...
%            C_gam_alp_beta{gamma,beta,alpha} ...
%        (I_g2 \otimes ZR(s3a_dum))
%
% INPUTS
%   G.C      : coefficient matrix CG, size g1*NG by g2
%   G.Z      : 1-by-ns3a cell array of exponent vectors
%   idxbeta  : nbeta-by-ns3a array, beta indices, entries in {1,2,3}
%   idxalpha : nalpha-by-ns3a array, alpha indices, entries in {1,2,3}
%   lims     : ns3a-by-2 domain array
%   Csize    : (optional) size of final C.params over gamma, used only for
%              an input consistency check
%
% OUTPUTS
%   C_gam_alp_beta : cell array of size 3^ns3a by nbeta by nalpha
%   ZL             : left monomial basis in s3a
%   ZR             : right monomial basis in s3a_dum
%
% Each C_gam_alp_beta{k,i,j} has size
%
%   g1*prod(cellfun(@numel,ZL)) by g2*prod(cellfun(@numel,ZR))
%
% The index convention for beta, alpha and gamma is
%   1 <-> 0   (multiplier, I_0(s) = delta(s))
%   2 <-> +1  (lower integral, I_1(s)  = 1 for s>=0)
%   3 <-> -1  (upper integral, I_-1(s) = 1 for s<=0)
%
% NOTES
% This is the shared version of the routine. Near-identical copies exist as
% @sopvar/private/int_semisep_AT.m and as local subfunctions of
% @sopvar/mtimes.m and @sopvar/mtimes_AT.m; those should eventually be
% replaced by calls to this function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int_semisep
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% AT, 2026: Initial coding as @sopvar/private/int_semisep_AT
% MP, 08/22/2026: Promoted to a shared function. Added a guard so that
%                 partial (nbeta<3^ns3a, nalpha<3^ns3a) index lists are
%                 supported, as the documented interface implies; the
%                 previous version passed a zero subscript to sub2ind
%                 whenever an enumerated key was absent from idxbeta or
%                 idxalpha.

ZG = G.Z(:).';
CG = G.C;

ns3a = numel(ZG);

% Normalize empty index arrays
if isempty(idxbeta)
    idxbeta = zeros(1,ns3a);
end
if isempty(idxalpha)
    idxalpha = zeros(1,ns3a);
end

nbeta  = size(idxbeta,1);
nalpha = size(idxalpha,1);

% Basic checks
if size(idxbeta,2) ~= ns3a
    error('int_semisep: idxbeta must have one column per s3a variable.');
end

if size(idxalpha,2) ~= ns3a
    error('int_semisep: idxalpha must have one column per s3a variable.');
end

if size(lims,1) ~= ns3a || size(lims,2) ~= 2
    error('int_semisep: lims must be numel(G.Z)-by-2.');
end

if any(idxbeta(:) < 1) || any(idxbeta(:) > 3) || ...
        any(idxbeta(:) ~= round(idxbeta(:)))
    error('int_semisep: idxbeta entries must be 1, 2, or 3.');
end

if any(idxalpha(:) < 1) || any(idxalpha(:) > 3) || ...
        any(idxalpha(:) ~= round(idxalpha(:)))
    error('int_semisep: idxalpha entries must be 1, 2, or 3.');
end

NG = prod(cellfun(@numel,ZG));

if NG == 0
    NG = 1;
end

if mod(size(CG,1),NG) ~= 0
    error('int_semisep: size(G.C,1) must be divisible by prod(numel(G.Z)).');
end

g1 = size(CG,1)/NG;
g2 = size(CG,2);

% No common pass-through variables.
% Then there is no semiseparable integration to perform.
if ns3a == 0
    C_gam_alp_beta = cell(1,nbeta,nalpha);

    for i = 1:nbeta
        for j = 1:nalpha
            C_gam_alp_beta{1,i,j} = CG;
        end
    end

    ZL = {};
    ZR = {};
    return
end

if nargin >= 5 && ~isempty(Csize) && prod(Csize) ~= 3^ns3a
    error('int_semisep: Csize is inconsistent with the number of s3a variables.');
end

a = lims(:,1);
b = lims(:,2);

% Build common left/right bases.
% These must support:
%   G(s), G(t), int G from a to s/t, int G from s/t to b,
%   and int G between s and t.
ZL = cell(1,ns3a);

for i = 1:ns3a
    E = ZG{i}(:);

    if any(E == -1)
        error('int_semisep: exponent -1 cannot be integrated by polynomial monomial rules.');
    end

    ZL{i} = unique([0; E; E+1]);
end

ZR = ZL;

NL = prod(cellfun(@numel,ZL));
NR = NL;

% Gamma multi-indices in the same linear order as 3-by-3-by-... parameter cells.
gamIdx = fliplr(dec2base(0:3^ns3a-1,3,ns3a)-'0') + 1;

C_gam_alp_beta = cell(3^ns3a,nbeta,nalpha);
for idx_C = 1:numel(C_gam_alp_beta)
    C_gam_alp_beta{idx_C} = sparse(g1*NL,g2*NR);
end

% The eleven (gamma,beta,alpha) keys for which p(gamma,alpha,beta) is
% nonzero, i.e. the eleven cases of the integral table in Sec. 6.1.
cllA = [111   212   313   221   331   222   232   332   223   323   333];

Ci_key_ell = cell(numel(cllA),ns3a);

for key_idx = 1:length(cllA)
    key = cllA(key_idx);

    for ell = 1:ns3a
        E = ZG{ell}(:);
        ZE = ZL{ell}(:);

        nE = numel(E);
        nZ = numel(ZE);

        Ci = sparse(nE,nZ^2);
        e = E;
        switch key

            % Evaluation at s
            case {111,212,313}
                coeffs = 0*e + 1;
                exp_s  = e;
                exp_t  = 0*e;

            % Evaluation at t
            case {221,331}
                coeffs = 0*e + 1;
                exp_s  = 0*e;
                exp_t  = e;

            % int_t^s eta^e deta
            case 222
                coeffs = [1./(e+1), -1./(e+1)];
                exp_s  = [e+1, 0*e];
                exp_t  = [0*e,   e+1];

            % int_s^b eta^e deta
            case 232
                coeffs = [b(ell).^(e+1)./(e+1), -1./(e+1)];
                exp_s  = [0*e, e+1];
                exp_t  = [0*e, 0*e];

            % int_t^b eta^e deta
            case 332
                coeffs = [b(ell).^(e+1)./(e+1), -1./(e+1)];
                exp_s  = [0*e, 0*e];
                exp_t  = [0*e, e+1];

            % int_a^t eta^e deta
            case 223
                coeffs = [-a(ell).^(e+1)./(e+1), 1./(e+1)];
                exp_s  = [0*e, 0*e];
                exp_t  = [0*e, e+1];

            % int_a^s eta^e deta
            case 323
                coeffs = [-a(ell).^(e+1)./(e+1), 1./(e+1)];
                exp_s  = [0*e, e+1];
                exp_t  = [0*e, 0*e];

            % int_s^t eta^e deta
            case 333
                coeffs = [-1./(e+1), 1./(e+1)];
                exp_s  = [e+1, 0*e];
                exp_t  = [0*e,   e+1];

            % Empty interval / incompatible ordering
            otherwise
                coeffs = [];
                exp_s  = [];
                exp_t  = [];
        end
        for h = 1:size(coeffs, 2)
            [is, ~, ~] = find(ZE == exp_s(:, h)');
            [it, ~, ~] = find(ZE == exp_t(:, h)');

            if isempty(is) || isempty(it)
                error('int_semisep: internal basis construction error.');
            end

            col = it + (is-1)*nZ;
            ind_Ci = sub2ind(size(Ci),1:length(e), col')';
            Ci(ind_Ci) = coeffs(:, h);
        end

        Ci_key_ell{key_idx, ell} = Ci;
    end
end

cllA = str2double(num2cell(num2str(cllA')));
all_possible_keys = cllA;

for ell = 1:(ns3a-1)
    all_possible_keys = [ kron(all_possible_keys, ones(size(cllA, 1), 1)), kron(ones(size(all_possible_keys, 1), 1), cllA)];
end

for key_idx_arg = 1:size(all_possible_keys, 1)
    chosen_gamma = all_possible_keys(key_idx_arg, 1:3:size(all_possible_keys, 2));
    chosen_beta  = all_possible_keys(key_idx_arg, 2:3:size(all_possible_keys, 2));
    chosen_alpha = all_possible_keys(key_idx_arg, 3:3:size(all_possible_keys, 2));

    [~, un_indices_alpha] = ismember(chosen_alpha, idxalpha, 'rows'); % j indices in C_gam_alp_beta
    [~, un_indices_beta]  = ismember(chosen_beta, idxbeta, 'rows');   % i indices in C_gam_alp_beta
    [~, un_indices_gamma] = ismember(chosen_gamma, gamIdx, 'rows');   % k indices in C_gam_alp_beta

    % This (gamma,beta,alpha) key is not requested by the caller.
    if un_indices_alpha==0 || un_indices_beta==0 || un_indices_gamma==0
        continue
    end

    Csep = 1;
    keys = all_possible_keys(key_idx_arg, :);
    for ell = 1:ns3a
        lfi = ((ell - 1) * 3 + 1);
        rgi = (ell*3 );
        key = keys(lfi:rgi);

        [~, key_idx] = ismember(key, cllA, 'rows');
        if key_idx==0
            error('int_semisep: internal key lookup error.');
        end
        Ci = Ci_key_ell{key_idx, ell};
        Csep = kron(Csep,Ci);
    end

    if ~isequal(size(Csep),[NG,NL*NR])
        error('int_semisep: internal Csep size mismatch.');
    end

    % Convert
    %
    %   (I_g1 \otimes Zmix')*(I_g1 \otimes Csep')*CG
    %
    % into
    %
    %   (I_g1 \otimes ZL')*Cnew*(I_g2 \otimes ZR).
    Cmul = kron(speye(g1),Csep.');
    Cnew = rearrangeCoef(Cmul,CG,g1,g2,ZL);
    if ~isequal(size(Cnew),[g1*NL,g2*NR])
        error('int_semisep: output coefficient size mismatch.');
    end

    ind_Cgam_alpha = sub2ind(size(C_gam_alp_beta), un_indices_gamma, un_indices_beta, un_indices_alpha);
    C_gam_alp_beta{ind_Cgam_alpha} = Cnew;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cout = rearrangeCoef(Cmul,C,p,q,Zint)
% rearrangeCoef
%
% Sparse-safe version.
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
% OUTPUT
%   Cout : sparse matrix of size p*NZ by q*NZ

Zint = Zint(:).';
ns = cellfun(@numel,Zint);

% Force sparse multiplication/storage.
D = sparse(Cmul) * sparse(C);

if isempty(ns)
    Cout = D;
    return
end

NZ = prod(ns);

if ~isequal(size(D),[p*NZ^2,q])
    error('rearrangeCoef: input dimensions are inconsistent.');
end

n = numel(ns);

% D rows are ordered as:
%
%   t_n, s_n, t_{n-1}, s_{n-1}, ..., t_1, s_1, p
%
% We want Cout ordered as:
%
%   rows: s_n, ..., s_1, p
%   cols: t_n, ..., t_1, q
%
% Instead of full(D), reshape, permute, reshape, we remap nonzero indices.

[rowD,colD,valD] = find(D);

oldDims = [repelem(fliplr(ns),2), p];

oldSub = cell(1,2*n+1);
[oldSub{:}] = ind2sub(oldDims,rowD);

% Extract old indices.
% oldSub = {t_n, s_n, t_{n-1}, s_{n-1}, ..., t_1, s_1, p}
tSub = oldSub(1:2:2*n);
sSub = oldSub(2:2:2*n);
pSub = oldSub{2*n+1};

% New row index: s_n, ..., s_1, p
rowDims = [fliplr(ns), p];
rowArgs = [sSub, {pSub}];
rowCout = sub2ind(rowDims,rowArgs{:});

% New column index: t_n, ..., t_1, q
colDims = [fliplr(ns), q];
colArgs = [tSub, {colD}];
colCout = sub2ind(colDims,colArgs{:});

Cout = sparse(rowCout,colCout,valD,NZ*p,NZ*q);

end

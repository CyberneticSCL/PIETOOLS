function [C_gam_alp_beta,ZL,ZR] = int_semisep(G,idxbeta,idxalpha,lims,Csize,layout)
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
%              an input consistency check. Pass [] to skip it
%   layout   : (optional) 'separate' (default) or 'packed', selecting the
%              output layout described below. May be given in place of
%              Csize, since the two are told apart by type
%
% OUTPUTS
%   C_gam_alp_beta : coefficient matrices, in one of two layouts
%   ZL             : left monomial basis in s3a
%   ZR             : right monomial basis in s3a_dum
%
% LAYOUTS
% 'separate' (default) keeps the beta and alpha axes as cell dimensions:
%
%   C_gam_alp_beta : 3^ns3a by nbeta by nalpha cell array, each entry
%                    g1*NL by g2*NR
%
% 'packed' collapses them into one matrix per gamma:
%
%   C_gam_alp_beta : 3^ns3a by 1 cell array, each entry
%                    g1*NL*nbeta by g2*NR*nalpha, holding block (i,j) at
%                    rows (1:g1*NL)+g1*NL*(i-1) and columns
%                    (1:g2*NR)+g2*NR*(j-1)
%
% where NL = prod(cellfun(@numel,ZL)) and NR = prod(cellfun(@numel,ZR)).
% Both hold the same numbers. 'packed' allocates 3^ns3a sparse matrices
% instead of 3^ns3a*nbeta*nalpha of them, which is cheaper when nbeta and
% nalpha are large, and is the form '@sopvar/mtimes_AT' consumes.
%
% The index convention for beta, alpha and gamma is
%   1 <-> 0   (multiplier, I_0(s) = delta(s))
%   2 <-> +1  (lower integral, I_1(s)  = 1 for s>=0)
%   3 <-> -1  (upper integral, I_-1(s) = 1 for s<=0)
%
% NOTES
% This is the only implementation of the routine. It is called directly by  % MMP, 08/30/2026
% 'possopvar' and '@sopvar/mtimes', and through the repacking wrapper       % MMP, 08/30/2026
% '@sopvar/private/int_semisep_AT' by '@sopvar/mtimes_AT'. It previously    % MMP, 08/30/2026
% existed in four near-identical copies, so a fix applied to one silently   % MMP, 08/30/2026
% missed the others.                                                        % MMP, 08/30/2026
%
% A repeated row in IDXBETA or IDXALPHA is filled in every slot it occupies,
% not only the first. Both callers pass deduplicated indices, so this costs
% nothing in practice, but the output does not depend on their doing so.

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

% % % Resolve the optional arguments. The layout may be supplied in either    % MMP, 08/30/2026
% % % trailing position, since Csize is numeric and layout is text.           % MMP, 08/30/2026
if nargin<5,    Csize = [];     end                                          % MMP, 08/30/2026
if nargin<6,    layout = '';    end                                          % MMP, 08/30/2026
if ischar(Csize) || isstring(Csize)                                          % MMP, 08/30/2026
    layout = Csize;     Csize = [];                                          % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026
if isempty(layout)                                                           % MMP, 08/30/2026
    layout = 'separate';                                                     % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026
if ~(ischar(layout) || isstring(layout)) ...                                 % MMP, 08/30/2026
        || ~any(strcmpi(layout,{'separate','packed'}))                       % MMP, 08/30/2026
    error('int_semisep: layout must be ''separate'' or ''packed''.');        % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026
packed = strcmpi(layout,'packed');                                           % MMP, 08/30/2026

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
    if packed                                                                % MMP, 08/30/2026
        C_gam_alp_beta = {repmat(CG,nbeta,nalpha)};                          % MMP, 08/30/2026
    else                                                                     % MMP, 08/30/2026
        C_gam_alp_beta = cell(1,nbeta,nalpha);
        for i = 1:nbeta
            for j = 1:nalpha
                C_gam_alp_beta{1,i,j} = CG;
            end
        end
    end                                                                      % MMP, 08/30/2026

    ZL = {};
    ZR = {};
    return
end

if ~isempty(Csize) && prod(Csize) ~= 3^ns3a                                  % MMP, 08/30/2026
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

% % % Precompute the reindexing that 'rearrangeCoef' applies.                 % MMP, 08/30/2026
%                                                                            % MMP, 08/30/2026
% Each key produces a matrix D whose row index encodes the multi-index        % MMP, 08/30/2026
% (t_n,s_n,...,t_1,s_1,p), and the result is that same nonzero moved to row   % MMP, 08/30/2026
% (s_n,...,s_1,p) and column (t_n,...,t_1,q). That permutation depends only   % MMP, 08/30/2026
% on the monomial sizes and on g1, not on the key, so computing it once here  % MMP, 08/30/2026
% replaces an ind2sub with 2*ns3a+1 outputs and two sub2ind calls per key by  % MMP, 08/30/2026
% two gathers. The scatter, not the multiply, is 69 to 89 percent of          % MMP, 08/30/2026
% 'rearrangeCoef'.                                                            % MMP, 08/30/2026
nsL = cellfun(@numel,ZL);                                                    % MMP, 08/30/2026
oldDims = [repelem(fliplr(nsL),2), g1];                                      % MMP, 08/30/2026
nrowD = g1*NL*NR;                                                            % MMP, 08/30/2026
oldSub = cell(1,2*ns3a+1);                                                   % MMP, 08/30/2026
[oldSub{:}] = ind2sub(oldDims,(1:nrowD)');                                   % MMP, 08/30/2026
tSub = oldSub(1:2:2*ns3a);                                                   % MMP, 08/30/2026
sSub = oldSub(2:2:2*ns3a);                                                   % MMP, 08/30/2026
rowMap = sub2ind([fliplr(nsL), g1], sSub{:}, oldSub{2*ns3a+1});              % MMP, 08/30/2026
if ns3a==1                                                                   % MMP, 08/30/2026
    colMap = tSub{1};                                                        % MMP, 08/30/2026
else                                                                         % MMP, 08/30/2026
    colMap = sub2ind(fliplr(nsL), tSub{:});                                  % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026

% Gamma multi-indices in the same linear order as 3-by-3-by-... parameter cells.
gamIdx = fliplr(dec2base(0:3^ns3a-1,3,ns3a)-'0') + 1;

% 'packed' needs 3^ns3a matrices rather than 3^ns3a*nbeta*nalpha of them.     % MMP, 08/30/2026
if packed                                                                    % MMP, 08/30/2026
    C_gam_alp_beta = cell(3^ns3a,1);                                         % MMP, 08/30/2026
    for idx_C = 1:numel(C_gam_alp_beta)                                      % MMP, 08/30/2026
        C_gam_alp_beta{idx_C} = sparse(g1*NL*nbeta,g2*NR*nalpha);            % MMP, 08/30/2026
    end                                                                      % MMP, 08/30/2026
else                                                                         % MMP, 08/30/2026
    C_gam_alp_beta = cell(3^ns3a,nbeta,nalpha);
    for idx_C = 1:numel(C_gam_alp_beta)
        C_gam_alp_beta{idx_C} = sparse(g1*NL,g2*NR);
    end
end                                                                          % MMP, 08/30/2026

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

% % % Enumerate only the (gamma,beta,alpha) triples the caller asked for.     % MMP, 08/30/2026
%                                                                            % MMP, 08/30/2026
% Per spatial direction the table holds 11 triples, so enumerating all of     % MMP, 08/30/2026
% them costs 11^ns3a and, when IDXBETA and IDXALPHA are the full 3^ns3a       % MMP, 08/30/2026
% index sets, that is exactly the number of nonzero output blocks and so is   % MMP, 08/30/2026
% optimal. It is badly pessimistic otherwise: a caller asking for one beta    % MMP, 08/30/2026
% row and one alpha row -- which is how 'possopvar' calls this routine --     % MMP, 08/30/2026
% needs at most 2^ns3a of them, since only (beta,alpha) = (3,2) and (2,3)     % MMP, 08/30/2026
% admit two gammas and every other pair admits one. Enumerating 11^ns3a and   % MMP, 08/30/2026
% discarding the rest therefore wasted a factor of 166 at ns3a = 3 and 5033   % MMP, 08/30/2026
% at ns3a = 5, which matters because this class is meant to take an           % MMP, 08/30/2026
% arbitrary number of spatial variables.                                      % MMP, 08/30/2026
%                                                                            % MMP, 08/30/2026
% Driving the loop from the requested rows instead makes the cost             % MMP, 08/30/2026
% proportional to the output actually asked for, and removes the need to      % MMP, 08/30/2026
% look beta and alpha back up: they are the loop indices. It also fills       % MMP, 08/30/2026
% repeated index rows for free.                                               % MMP, 08/30/2026
pw3 = 3.^(0:ns3a-1)';                                                        % MMP, 08/30/2026
key_of = cell(3,3);                                                          % MMP, 08/30/2026
for k = 1:size(cllA,1)                                                       % MMP, 08/30/2026
    key_of{cllA(k,2),cllA(k,3)} = [key_of{cllA(k,2),cllA(k,3)}, k];          % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026

for ib = 1:nbeta                                                             % MMP, 08/30/2026
for ia = 1:nalpha                                                            % MMP, 08/30/2026
    % Options per direction, and the size of their Cartesian product.        % MMP, 08/30/2026
    opt = cell(1,ns3a);     nopt = zeros(1,ns3a);                            % MMP, 08/30/2026
    for ell = 1:ns3a                                                         % MMP, 08/30/2026
        opt{ell} = key_of{idxbeta(ib,ell),idxalpha(ia,ell)};                  % MMP, 08/30/2026
        nopt(ell) = numel(opt{ell});                                         % MMP, 08/30/2026
    end                                                                      % MMP, 08/30/2026
    if any(nopt==0)                                                          % MMP, 08/30/2026
        continue                                                             % MMP, 08/30/2026
    end                                                                      % MMP, 08/30/2026

for comb = 0:prod(nopt)-1                                                    % MMP, 08/30/2026
    un_indices_beta  = ib;                                                   % MMP, 08/30/2026
    un_indices_alpha = ia;                                                   % MMP, 08/30/2026

    % Mixed-radix walk over the Cartesian product of the per-direction       % MMP, 08/30/2026
    % options, accumulating both the Kronecker factor and the linear index    % MMP, 08/30/2026
    % of gamma. The gamma multi-index encodes to its own cell position in     % MMP, 08/30/2026
    % base 3, so no lookup table is needed.                                   % MMP, 08/30/2026
    Csep = 1;   gam_lin = 1;    rem_c = comb;                                % MMP, 08/30/2026
    for ell = 1:ns3a                                                         % MMP, 08/30/2026
        sel = mod(rem_c,nopt(ell)) + 1;                                      % MMP, 08/30/2026
        rem_c = floor(rem_c/nopt(ell));                                      % MMP, 08/30/2026
        key_idx = opt{ell}(sel);                                             % MMP, 08/30/2026
        Csep = kron(Csep,Ci_key_ell{key_idx,ell});                           % MMP, 08/30/2026
        gam_lin = gam_lin + (cllA(key_idx,1)-1)*pw3(ell);                    % MMP, 08/30/2026
    end                                                                      % MMP, 08/30/2026
    un_indices_gamma = gam_lin;                                              % MMP, 08/30/2026

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
    Cnew = rearrangeCoef(Csep,CG,g1,g2,rowMap,colMap,NL);                   % MMP, 08/30/2026
    if ~isequal(size(Cnew),[g1*NL,g2*NR])
        error('int_semisep: output coefficient size mismatch.');
    end

    if packed                                                                % MMP, 08/30/2026
        rows = (1:(g1*NL)) + g1*NL*(un_indices_beta - 1);                    % MMP, 08/30/2026
        cols = (1:(g2*NR)) + g2*NR*(un_indices_alpha - 1);                   % MMP, 08/30/2026
        C_gam_alp_beta{un_indices_gamma}(rows,cols) = Cnew;                  % MMP, 08/30/2026
    else                                                                     % MMP, 08/30/2026
        ind_Cgam_alpha = sub2ind(size(C_gam_alp_beta), un_indices_gamma, un_indices_beta, un_indices_alpha);
        C_gam_alp_beta{ind_Cgam_alpha} = Cnew;
    end                                                                      % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026

% Repeated rows in IDXBETA or IDXALPHA need no special handling: the loop     % MMP, 08/30/2026
% above is driven by the row indices themselves, so every occurrence is       % MMP, 08/30/2026
% filled. The earlier version resolved a multi-index back to its first        % MMP, 08/30/2026
% matching row and needed a separate pass to copy the duplicates.             % MMP, 08/30/2026

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cout = rearrangeCoef(Csep,C,p,q,rowMap,colMap,NZ)                   % MMP, 08/30/2026
% rearrangeCoef
%
% Converts
%
%   (I_p o Zmix(s,t)') * kron(I_p,Csep') * C
%
% into
%
%   (I_p o ZL(s)') * Cout * (I_q o ZR(t))
%
% by moving each stored nonzero to the position the second form wants. Both
% steps are cheap:
%
%   D = kron(I_p,Csep')*C is evaluated blockwise. kron(I_p,M) is block
%   diagonal, so materializing it would cost p times the nonzeros of M and
%   make the multiply p times wider than it need be; reshaping C so that its
%   p row blocks sit side by side gives the same result from one multiply.
%
%   The reindexing is a fixed permutation of D's row index, depending only on
%   the monomial sizes and on p, not on the key being processed. ROWMAP and
%   COLMAP carry it, precomputed once by the caller, so what remains here is
%   two gathers rather than an ind2sub with 2*ns3a+1 outputs and two sub2ind
%   calls per key. The scatter, not the multiply, was 69 to 89 percent of
%   this routine.
%
% INPUTS
%   Csep   : NG by NZ^2 separable coefficient matrix for one key
%   C      : p*NG by q coefficient matrix of G
%   p,q    : row and column block counts
%   rowMap : p*NZ^2 vector, destination row for each row of D
%   colMap : p*NZ^2 vector, destination column of D's row index, before the
%            contribution of D's own column index
%   NZ     : prod(cellfun(@numel,ZL))
%
% OUTPUT
%   Cout   : NZ*p by NZ*q sparse matrix

D = reshape(sparse(Csep).' * reshape(sparse(C),size(Csep,1),[]), [], q);

if isempty(rowMap)
    Cout = D;
    return
end
if size(D,1)~=numel(rowMap)
    error('rearrangeCoef: input dimensions are inconsistent.');
end

[rowD,colD,valD] = find(D);
Cout = sparse(rowMap(rowD), colMap(rowD) + NZ*(colD-1), valD, NZ*p, NZ*q);

end

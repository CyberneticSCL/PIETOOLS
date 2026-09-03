% RETIRED 2026-08-30. Local subfunction formerly in @sopvar/mtimes.m
% (lines 666-1015), removed when mtimes was pointed at the shared
% sopvar/int_semisep.m. Verified equivalent to that routine except that it
% fills every repeated row of idxbeta/idxalpha rather than only the first;
% the shared routine now does the same. Kept for reference only: this
% folder is inside an @-folder and so is never on the MATLAB path.

function [C_gam_alp_beta,ZL,ZR] = int_semisep(G,idxbeta,idxalpha,lims,Csize)
% int_semisep
%
% Converts
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
%   idxbeta  : nbeta-by-ns3a array, beta_b indices
%   idxalpha : nalpha-by-ns3a array, alpha_a indices
%   lims     : ns3a-by-2 domain array
%   Csize    : size of final C.params over gamma, used for compatibility
%
% OUTPUTS
%   C_gam_alp_beta : cell array of size 3^ns3a by nbeta by nalpha
%   ZL             : left monomial basis in s3a
%   ZR             : right monomial basis in s3a_dum
%
% Each C_gam_alp_beta{k,i,j} has size
%
%   g1*prod(cellfun(@numel,ZL)) by g2*prod(cellfun(@numel,ZR))

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

if nargin >= 5 && prod(Csize) ~= 3^ns3a
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

% Main loop over beta, alpha, gamma
for i = 1:nbeta

    beta = idxbeta(i,:);

    for j = 1:nalpha

        alpha = idxalpha(j,:);

        for k = 1:3^ns3a

            gamma = gamIdx(k,:);

            % Csep maps original ZG(eta) into the separated mixed basis
            % kron_i (ZL_i(s_i) \otimes ZR_i(t_i)).
            %
            % Csep has size
            %
            %   NG by NL*NR
            Csep = 1;

            for ell = 1:ns3a

                E = ZG{ell}(:);
                ZE = ZL{ell}(:);

                nE = numel(E);
                nZ = numel(ZE);

                % Local basis is kron(ZL_ell(s),ZR_ell(t)).
                % Column index is:
                %
                %   col = idx_t + (idx_s-1)*nZ
                %
                % so Ci has size nE by nZ^2.
                Ci = sparse(nE,nZ^2);

                key = 100*gamma(ell) + 10*beta(ell) + alpha(ell);

                for r = 1:nE

                    e = E(r);

                    switch key

                        % Evaluation at s
                        case {111,212,313}
                            coeffs = 1;
                            exp_s  = e;
                            exp_t  = 0;

                        % Evaluation at t
                        case {221,331}
                            coeffs = 1;
                            exp_s  = 0;
                            exp_t  = e;

                        % int_t^s eta^e deta
                        case 222
                            coeffs = [1/(e+1), -1/(e+1)];
                            exp_s  = [e+1, 0];
                            exp_t  = [0,   e+1];

                        % int_s^b eta^e deta
                        case 232
                            coeffs = [b(ell)^(e+1)/(e+1), -1/(e+1)];
                            exp_s  = [0, e+1];
                            exp_t  = [0, 0];

                        % int_t^b eta^e deta
                        case 332
                            coeffs = [b(ell)^(e+1)/(e+1), -1/(e+1)];
                            exp_s  = [0, 0];
                            exp_t  = [0, e+1];

                        % int_a^t eta^e deta
                        case 223
                            coeffs = [-a(ell)^(e+1)/(e+1), 1/(e+1)];
                            exp_s  = [0, 0];
                            exp_t  = [0, e+1];

                        % int_a^s eta^e deta
                        case 323
                            coeffs = [-a(ell)^(e+1)/(e+1), 1/(e+1)];
                            exp_s  = [0, e+1];
                            exp_t  = [0, 0];

                        % int_s^t eta^e deta
                        case 333
                            coeffs = [-1/(e+1), 1/(e+1)];
                            exp_s  = [e+1, 0];
                            exp_t  = [0,   e+1];

                        % Empty interval / incompatible ordering
                        otherwise
                            coeffs = [];
                            exp_s  = [];
                            exp_t  = [];
                    end

                    for h = 1:numel(coeffs)
                        is = find(ZE == exp_s(h),1);
                        it = find(ZE == exp_t(h),1);

                        if isempty(is) || isempty(it)
                            error('int_semisep: internal basis construction error.');
                        end

                        col = it + (is-1)*nZ;
                        Ci(r,col) = Ci(r,col) + coeffs(h);
                    end
                end

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

            C_gam_alp_beta{k,i,j} = Cnew;
        end
    end
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

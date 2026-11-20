function [C_gam_alp_beta,ZL,ZR] = int_semisep(G,idxbeta,idxalpha,lims,Csize,layout)
% [C_GAM_ALP_BETA,ZL,ZR] = INT_SEMISEP(G,IDXBETA,IDXALPHA,LIMS,CSIZE,LAYOUT)
% evaluates the two-sided semiseparable integral
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
%   idxbeta  : nbeta-by-ns3a array, beta indices, entries in {1,2,3}, any
%              subset of the multi-indices
%   idxalpha : nalpha-by-ns3a array, alpha indices, entries in {1,2,3},
%              any subset of the multi-indices
%   lims     : ns3a-by-2 domain array
%   Csize    : (optional) size of final C.params over gamma, used only for
%              an input consistency check. Pass [] to skip it
%   layout   : (optional) 'separate' (default), 'packed' or 'triplets',
%              selecting the output layout described below. May be given
%              in place of Csize, since the two are told apart by type
%
% OUTPUTS
%   C_gam_alp_beta : coefficient matrices, in one of three layouts
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
% nalpha are large.
%
% 'triplets' is 'separate' with each matrix M replaced by the struct
%   struct('i',I,'j',J,'v',V,'m',size(M,1),'n',size(M,2))
% of column vectors, M == sparse(I,J,V,m,n), no repeated (I,J), no zero V,
% in no particular order. The matrix is never formed: its g2*NR columns can
% carry the decision variables, so its column pointers alone would cost
% O(q) per cell to write and again to scan.
%
% The index convention for beta, alpha and gamma is
%   1 <-> 0   (multiplier, I_0(s) = delta(s))
%   2 <-> +1  (lower integral, I_1(s)  = 1 for s>=0)
%   3 <-> -1  (upper integral, I_-1(s) = 1 for s<=0)
%
% NOTES
% A repeated row in IDXBETA or IDXALPHA is filled in every slot it occupies,
% not only the first.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

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
% AT, 2026: Initial coding

% % % Resolve the optional arguments. The layout may be supplied in either
% % % trailing position, since Csize is numeric and layout is text.
if nargin<5,    Csize = [];     end
if nargin<6,    layout = '';    end
if ischar(Csize) || isstring(Csize)
    layout = Csize;     Csize = [];
end
if isempty(layout)
    layout = 'separate';
end
if ~(ischar(layout) || isstring(layout)) ...
        || ~any(strcmpi(layout,{'separate','packed','triplets'}))
    error(['int_semisep: layout must be ''separate'', ''packed'' or ' ...
           '''triplets''.']);
end
packed = strcmpi(layout,'packed');
trip   = strcmpi(layout,'triplets');    % 'separate' shape, triplet entries

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

% One domain row per variable; index entries 1, 2 or 3
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

% Basis size of G, and the row and column block counts of G.C
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
    if packed
        C_gam_alp_beta = {repmat(CG,nbeta,nalpha)};
    else
        if trip     % CG itself is the output; hand over its nonzeros
            [ti,tj,tv] = find(CG);
            CG = struct('i',ti(:),'j',tj(:),'v',tv(:), ...
                        'm',size(CG,1),'n',size(CG,2));
        end
        % CG in every (beta,alpha) cell
        C_gam_alp_beta = cell(1,nbeta,nalpha);
        for i = 1:nbeta
            for j = 1:nalpha
                C_gam_alp_beta{1,i,j} = CG;
            end
        end
    end

    ZL = {};
    ZR = {};
    return
end

if ~isempty(Csize) && prod(Csize) ~= 3^ns3a
    error('int_semisep: Csize is inconsistent with the number of s3a variables.');
end

% Three steps:
%   setup     the bases ZL, ZR, the reindexing of 'rearrangeCoef' and the
%             per-direction factors of the eleven keys; it depends on
%             (G.Z, lims, g1) alone, so it is fetched from a bounded memo;
%   enumerate for each requested (beta,alpha) row pair, the gammas with
%             p(gamma,alpha,beta) nonzero and their factors;
%   assemble  the coefficients of each block, in the layout.
S = memo_setup(ZG,lims,g1);
ZL = S.ZL;      ZR = S.ZR;
[blk,kseqs] = semisep_enumerate(S,idxbeta,idxalpha);
C_gam_alp_beta = semisep_assemble(S,blk,kseqs,CG,NG,g1,g2,nbeta,nalpha, ...
                                  packed,trip);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = memo_setup(ZG,lims,g1)
% S = MEMO_SETUP(ZG,LIMS,G1) returns the setup of (ZG, lims, g1), from the
% memo when it holds that key. Calls with G.Z and the domain fixed come
% with a few distinct keys interleaved (measured on program builds: 9 or
% 10 keys, 54 to 68 key changes in up to 3251 calls), so the memo holds
% several setups, not the last one only. Keys match bit for bit; the last
% match is tried first. Bounded: at most 16 setups and 2^20 rowMap entries
% in all (rowMap and colMap have g1*NL*NR entries each: 16 MB), oldest
% dropped first. A setup is of spatial size and holds no G.C, so nothing
% q-sized is kept.
% INPUT   ZG: 1 x ns3a cell of exponent columns; lims: ns3a x 2; g1: the
%         row block count of G.C.
% OUTPUT  S: the setup struct (fields as returned by semisep_setup).
persistent keys Ss last
% The setup key (ZG, lims, g1) as the 64-bit patterns of g1, the counts and
% every entry, so equal keys give bit-identical setups; isequal on the
% values would equate 0 and -0. [] (never cached) unless lims and every
% ZG{i} are real full double and each ZG{i} a column: those fix the class
% and shape of the setup, and the counts then fix every size. lims is
% numel(ZG) x 2 (checked in semisep_setup). A few builtins, independent of
% numel(ZG): an entry-by-entry comparison cost 6-22 us per call (measured).
key = [];
if isa(lims,'double') && ~issparse(lims) && isreal(lims) ...
        && all(cellfun('isclass',ZG,'double')) && all(cellfun('isreal',ZG)) ...
        && all(cellfun('size',ZG,2)==1)
    v = vertcat(ZG{:});         % sparse if any ZG{i} is
    if ~issparse(v)
        key = typecast([g1; numel(ZG); cellfun('prodofsize',ZG(:)); ...
                        lims(:); v],'uint64');
    end
end
% The last match first, then every cached key
if ~isempty(key) && ~isempty(last) && isequal(key,keys{last})
    S = Ss{last};
    return
end
if ~isempty(key)
    for i = 1:numel(keys)
        if isequal(key,keys{i}),    S = Ss{i};  last = i;   return,     end
    end
end
% A miss: build the setup, and cache it within the bounds, oldest first out
S = semisep_setup(ZG,lims,g1);
if ~isempty(key) && numel(S.rowMap)<=2^20
    keys{end+1} = key;      Ss{end+1} = S;
    nmap = cellfun(@(s) numel(s.rowMap),Ss);
    while numel(Ss)>16 || sum(nmap)>2^20
        keys(1) = [];   Ss(1) = [];     nmap(1) = [];
    end
    last = numel(Ss);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = semisep_setup(ZG,lims,g1)
% S = SEMISEP_SETUP(ZG,LIMS,G1) returns everything of int_semisep that
% depends on the key (ZG, lims, g1) alone: the common bases ZL = ZR, the
% permutation of 'rearrangeCoef'
% (rowMap, colMap), the factor Ci_key_ell{key,ell} of each key and
% direction with the classes of equal factors (clsOf, nCls, pwCls), and
% the key lists key_of{beta,alpha}. Nothing here reads G.C, the index rows
% or the layout, which is what makes caching S exact.
% INPUT   ZG: 1 x ns3a cell of exponent columns; lims: ns3a x 2; g1: the
%         row block count of G.C.
% OUTPUT  S: struct with fields ZL, ZR, NL, NR, rowMap, colMap, Ci_key_ell,
%         cllA, clsOf, nCls, pwCls, pw3, key_of.
ns3a = numel(ZG);
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

% % % Precompute the reindexing that 'rearrangeCoef' applies: the row index
% of D encodes (t_n,s_n,...,t_1,s_1,p), and its nonzero moves to row
% (s_n,...,s_1,p) and column (t_n,...,t_1,q). The permutation depends on the
% monomial sizes and g1 only, not on the key, so it is two gathers per key.
nsL = cellfun(@numel,ZL);
oldDims = [repelem(fliplr(nsL),2), g1];
nrowD = g1*NL*NR;
oldSub = cell(1,2*ns3a+1);
[oldSub{:}] = ind2sub(oldDims,(1:nrowD)');
tSub = oldSub(1:2:2*ns3a);
sSub = oldSub(2:2:2*ns3a);
rowMap = sub2ind([fliplr(nsL), g1], sSub{:}, oldSub{2*ns3a+1});
if ns3a==1
    colMap = tSub{1};
else
    colMap = sub2ind(fliplr(nsL), tSub{:});
end

% The eleven (gamma,beta,alpha) keys for which p(gamma,alpha,beta) is
% nonzero, i.e. the eleven cases of the integral table.
cllA = [111   212   313   221   331   222   232   332   223   323   333];
% The same eleven keys as rows [gamma beta alpha], row k the digits of
% cllA(k): column 1 gamma, 2 beta, 3 alpha. cllA is set to it after the
% switch.
keyGBA = [1 1 1; 2 1 2; 3 1 3; 2 2 1; 3 3 1; 2 2 2; ...
          2 3 2; 3 3 2; 2 2 3; 3 2 3; 3 3 3];

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
        % Term h of each exponent e: coefficient coeffs(:,h) at the basis
        % position of s^exp_s * t^exp_t
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

% From here on, cllA holds the keys as rows [gamma beta alpha].
cllA = keyGBA;

% The eleven keys share only EIGHT distinct per-direction factors: the switch
% above handles {111,212,313} in one branch and {221,331} in another, so
% Ci_key_ell takes 8 values per direction and the Kronecker product Csep
% below takes 8^ns3a, not 11^ns3a. Every repeat costs a redundant
% 'rearrangeCoef', which is 11^ns3a/8^ns3a = 1.89x at two pass-through
% variables, 2.60x at three and 3.57x at four. The duplicates are found by
% COMPARING the built factors rather than by hard-coding the branch labels,
% so the grouping stays correct if the table or the domains change; the
% comparison is 11^2 per direction on matrices of a few dozen nonzeros.
nkey  = size(cllA,1);
repOf = repmat((1:nkey).',1,ns3a);
for ell = 1:ns3a
    for k = 2:nkey
        for j = 1:k-1
            if isequal(Ci_key_ell{k,ell},Ci_key_ell{j,ell})
                repOf(k,ell) = repOf(j,ell);
                break
            end
        end
    end
end
% Renumber the representatives densely per direction, so the memo is indexed
% in mixed radix over the DISTINCT counts (8 per direction in the current
% table) rather than over 11, and holds prod(nCls) slots instead of
% 11^ns3a.
clsOf = zeros(nkey,ns3a);   nCls = zeros(1,ns3a);
for ell = 1:ns3a
    [~,~,ic] = unique(repOf(:,ell));
    clsOf(:,ell) = ic;
    nCls(ell)    = max(ic);
end
pwCls = cumprod([1,nCls(1:end-1)]);

% Base-3 weights of gamma, and the keys of each (beta,alpha) pair
pw3 = 3.^(0:ns3a-1)';
key_of = cell(3,3);
for k = 1:size(cllA,1)
    key_of{cllA(k,2),cllA(k,3)} = [key_of{cllA(k,2),cllA(k,3)}, k];
end
S = struct('ZL',{ZL},'ZR',{ZR},'NL',NL,'NR',NR,'rowMap',rowMap, ...
           'colMap',colMap,'Ci_key_ell',{Ci_key_ell},'cllA',cllA, ...
           'clsOf',clsOf,'nCls',nCls,'pwCls',pwCls,'pw3',pw3, ...
           'key_of',{key_of});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [blk,kseqs] = semisep_enumerate(S,idxbeta,idxalpha)
% [BLK,KSEQS] = SEMISEP_ENUMERATE(S,IDXBETA,IDXALPHA) lists the blocks to
% build: p(gamma,alpha,beta) is the product over directions of the
% per-direction table, so for a (beta,alpha) row pair
% the gammas are the Cartesian product of the per-direction key options.
% Row t of blk is [gamma cell, beta row, alpha row, memo slot] of block t,
% in the order the loop below met it; kseqs(t,:) its key per direction.
% Index work only, on at most nbeta*nalpha*2^ns3a rows.
% INPUT   S: the setup struct; idxbeta, idxalpha: the index rows.
% OUTPUT  blk: nb x 4 block table; kseqs: nb x ns3a key indices.
ns3a = numel(S.ZL);
nbeta = size(idxbeta,1);        nalpha = size(idxalpha,1);
key_of = S.key_of;  cllA = S.cllA;  clsOf = S.clsOf;
pwCls = S.pwCls;    pw3 = S.pw3;
% A (beta,alpha) pair admits at most max numel(key_of) keys per direction.
nmax = nbeta*nalpha*max(cellfun('prodofsize',key_of(:)))^ns3a;
blk = zeros(nmax,4);    kseqs = zeros(nmax,ns3a);   nb = 0;

% % % Enumerate only the (gamma,beta,alpha) triples the caller asked for:
% a (beta,alpha) pair admits one gamma per direction, two for (3,2) and
% (2,3), so at most 2^ns3a blocks per requested row pair against the
% 11^ns3a of the full table (a factor 166 at ns3a = 3 for one row pair).
% The row indices are the loop indices, so repeated rows are filled as they
% come.
for ib = 1:nbeta
for ia = 1:nalpha
    % Options per direction, and the size of their Cartesian product.
    opt = cell(1,ns3a);     nopt = zeros(1,ns3a);
    for ell = 1:ns3a
        opt{ell} = key_of{idxbeta(ib,ell),idxalpha(ia,ell)};
        nopt(ell) = numel(opt{ell});
    end
    if any(nopt==0)
        continue
    end

for comb = 0:prod(nopt)-1
    un_indices_beta  = ib;
    un_indices_alpha = ia;

    % Mixed-radix walk over the Cartesian product of the per-direction
    % options, accumulating both the Kronecker factor and the linear index
    % of gamma. The gamma multi-index encodes to its own cell position in
    % base 3, so no lookup table is needed.
    % Walk the mixed radix for the keys, gamma, and the memo slot. Csep is
    % NOT built here: it is only needed on a memo miss, so the Kronecker
    % chain is skipped on a hit along with the rearrangeCoef.
    gam_lin = 1;    rem_c = comb;    memoLin = 1;
    kseq = zeros(1,ns3a);
    for ell = 1:ns3a
        sel = mod(rem_c,nopt(ell)) + 1;
        rem_c = floor(rem_c/nopt(ell));
        key_idx = opt{ell}(sel);
        kseq(ell) = key_idx;
        memoLin = memoLin + (clsOf(key_idx,ell)-1)*pwCls(ell);
        gam_lin = gam_lin + (cllA(key_idx,1)-1)*pw3(ell);
    end
    un_indices_gamma = gam_lin;
    % Record the block; semisep_assemble builds it.
    nb = nb + 1;
    blk(nb,:) = [un_indices_gamma,un_indices_beta, ...
                 un_indices_alpha,memoLin];
    kseqs(nb,:) = kseq;
end
end
end

% Repeated rows in IDXBETA or IDXALPHA need no special handling: the loop
% above is driven by the row indices themselves, so every occurrence is
% filled.
blk = blk(1:nb,:);      kseqs = kseqs(1:nb,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C_gam_alp_beta = semisep_assemble(S,blk,kseqs,CG,NG,g1,g2, ...
                                           nbeta,nalpha,packed,trip)
% C_GAM_ALP_BETA = SEMISEP_ASSEMBLE(S,BLK,KSEQS,CG,NG,G1,G2,NBETA,NALPHA,
% PACKED,TRIP) returns the coefficients of every block of blk: the
% Kronecker product of its per-direction factors kseqs, applied to CG and
% re-indexed onto ZL, ZR by 'rearrangeCoef', memoized per slot since blocks
% share factors; written into the requested layout. Blocks are taken in
% the order of blk, so each memo slot is built from its first block.
% INPUT   S, blk, kseqs: from the setup and the enumeration; CG: G.C; NG,
%         g1, g2, nbeta, nalpha: sizes; packed, trip: layout flags.
% OUTPUT  C_gam_alp_beta: the coefficients in the requested layout.
ns3a = numel(S.ZL);     NL = S.NL;      NR = S.NR;      nCls = S.nCls;
Ci_key_ell = S.Ci_key_ell;  rowMap = S.rowMap;  colMap = S.colMap;
% 'packed' needs 3^ns3a matrices rather than 3^ns3a*nbeta*nalpha of them.
if packed
    % The packed outputs are accumulated as triplets and built with one
    % 'sparse' call each, below the loop. One slot per
    % (gamma,beta,alpha): the eleven keys are
    % distinct triples and gamma -> gam_lin is a bijection, so each block
    % is written at most once and nothing needs summing.
    nGam = 3^ns3a;
    accI = cell(nGam,nbeta*nalpha);
    accJ = cell(nGam,nbeta*nalpha);
    accV = cell(nGam,nbeta*nalpha);
elseif trip
    % Unwritten cells: no triplets, the dimensions of the all-zero matrix.
    C_gam_alp_beta = cell(3^ns3a,nbeta,nalpha);
    C_gam_alp_beta(:) = {struct('i',zeros(0,1),'j',zeros(0,1), ...
        'v',zeros(0,1),'m',g1*NL,'n',g2*NR)};
else
    C_gam_alp_beta = cell(3^ns3a,nbeta,nalpha);
    for idx_C = 1:numel(C_gam_alp_beta)
        C_gam_alp_beta{idx_C} = sparse(g1*NL,g2*NR);
    end
end
memoC = cell(prod([nCls,1]),1);
for t = 1:size(blk,1)
    un_indices_gamma = blk(t,1);        un_indices_beta = blk(t,2);
    un_indices_alpha = blk(t,3);        memoLin = blk(t,4);
    kseq = kseqs(t,:);

    % Convert
    %
    %   (I_g1 \otimes Zmix')*(I_g1 \otimes Csep')*CG
    %
    % into
    %
    %   (I_g1 \otimes ZL')*Cnew*(I_g2 \otimes ZR).
    % Keyed on the representative tuple, so the 11^ns3a enumerated blocks
    % share the 8^ns3a distinct results. The size checks stay on the miss
    % path, so every distinct Csep and Cnew is still validated once.
    if isempty(memoC{memoLin})
        Csep = 1;
        for ell = 1:ns3a
            Csep = kron(Csep,Ci_key_ell{kseq(ell),ell});
        end
        if ~isequal(size(Csep),[NG,NL*NR])
            error('int_semisep: internal Csep size mismatch.');
        end
        % 'triplets': Cmiss is a triplet struct and its size is its m, n.
        Cmiss = rearrangeCoef(Csep,CG,g1,g2,rowMap,colMap,NL,trip);
        if trip,  szC = [Cmiss.m,Cmiss.n];  else,  szC = size(Cmiss);  end
        if ~isequal(szC,[g1*NL,g2*NR])
            error('int_semisep: output coefficient size mismatch.');
        end
        memoC{memoLin} = Cmiss;
    end
    Cnew = memoC{memoLin};

    if packed
        % Stash the nonzeros of this block, shifted onto its (beta,alpha)
        % offset.
        [ri,ci,vi] = find(Cnew);
        slot = un_indices_beta + nbeta*(un_indices_alpha - 1);
        accI{un_indices_gamma,slot} = ri + g1*NL*(un_indices_beta - 1);
        accJ{un_indices_gamma,slot} = ci + g2*NR*(un_indices_alpha - 1);
        accV{un_indices_gamma,slot} = vi;
    else
        ind_Cgam_alpha = sub2ind(size(C_gam_alp_beta), un_indices_gamma, un_indices_beta, un_indices_alpha);
        C_gam_alp_beta{ind_Cgam_alpha} = Cnew;
    end
end

% One 'sparse' per gamma from the accumulated triplets. A gamma that was
% never written concatenates to empty, which 'sparse' turns into the
% all-zero matrix.
if packed
    C_gam_alp_beta = cell(nGam,1);
    for k = 1:nGam
        C_gam_alp_beta{k} = sparse(vertcat(accI{k,:}),vertcat(accJ{k,:}), ...
            vertcat(accV{k,:}),g1*NL*nbeta,g2*NR*nalpha);
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cout = rearrangeCoef(Csep,C,p,q,rowMap,colMap,NZ,trip)
% COUT = REARRANGECOEF(CSEP,C,P,Q,ROWMAP,COLMAP,NZ,TRIP) converts
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
%   When that reshape is wide -- q carries the decision variables and its
%   p*q columns are nearly all empty -- only its nonempty columns are
%   formed, from the triplets of C.
%
%   The reindexing is a fixed permutation of the row index of D, depending
%   only on the monomial sizes and on p, not on the key being processed.
%   ROWMAP and COLMAP carry it, precomputed once per setup, so what remains
%   here is two gathers rather than an ind2sub with 2*ns3a+1 outputs and two
%   sub2ind calls per key.
%
% INPUTS
%   Csep   : NG by NZ^2 separable coefficient matrix for one key
%   C      : p*NG by q coefficient matrix of G
%   p,q    : row and column block counts
%   rowMap : p*NZ^2 vector, destination row for each row of D
%   colMap : p*NZ^2 vector, destination column for each row of D, before
%            the contribution of the column index of D
%   NZ     : prod(cellfun(@numel,ZL))
%   trip   : (optional, default false) return Cout as the triplet struct
%            of the 'triplets' layout instead of building it
%
% OUTPUT
%   Cout   : NZ*p by NZ*q sparse matrix
if nargin<8,    trip = false;   end

% Wide X: the general lines below form X = reshape(C,NG,[]) with p*q
% columns and the product over all of them, O(p*q) whatever nnz(C); when q
% carries the decision variables nearly every column is empty. Here only
% the nonempty columns of X are built, from the triplets of C, O(nnz(C)).
% Exact: sparse*sparse computes each column of the product from the same
% column of X alone, so the result is bit-identical. The crossover rule is
% measured (p = 8, 25, 64); a narrow X takes the general lines.
if ~isempty(rowMap) && p*q > 8*nnz(C) + 1024
    NG  = size(Csep,1);     NZ2 = size(Csep,2);
    if NZ2*size(C,1)~=NG*numel(rowMap)          % the size(D,1) test below
        error('rearrangeCoef: input dimensions are inconsistent.');
    end
    [ci,cj,cv] = find(C);
    ci = ci(:);     cj = cj(:);     cv = cv(:);     % rows if C is a row
    rr   = mod(ci-1,NG) + 1;                        % monomial index in G
    xcol = (ci-rr)/NG + 1 + p*(cj-1);               % column of X: pp+p*(c-1)
    % find runs column-major with rows ascending, so xcol is nondecreasing
    % and a run-length code numbers its distinct values without a sort.
    isnew = diff([0; xcol])~=0;
    ux    = xcol(isnew);
    Xc    = sparse(rr, cumsum(isnew), cv, NG, numel(ux));
    [kD,jD,valD] = find(sparse(Csep).' * Xc);
    pc    = ux(jD);                                 % column of X
    ppD   = mod(pc-1,p) + 1;
    rowD  = kD + NZ2*(ppD-1);                       % row of D
    colD  = (pc-ppD)/p + 1;                         % column of D
    if trip     % triplets of the scatter: bijective, so no duplicates
        Cout = struct('i',rowMap(rowD(:)),'j',colMap(rowD(:))+NZ*(colD(:)-1), ...
                      'v',valD(:),'m',NZ*p,'n',NZ*q);
        return
    end
    Cout = sparse(rowMap(rowD), colMap(rowD) + NZ*(colD-1), valD, NZ*p, NZ*q);
    return
end

% General case: D = kron(I_p,Csep')*C from one product, then the scatter
D = reshape(sparse(Csep).' * reshape(sparse(C),size(Csep,1),[]), [], q);

if isempty(rowMap)
    Cout = D;
    return
end
if size(D,1)~=numel(rowMap)
    error('rearrangeCoef: input dimensions are inconsistent.');
end

[rowD,colD,valD] = find(D);
if trip         % as in the wide branch above
    Cout = struct('i',rowMap(rowD(:)),'j',colMap(rowD(:))+NZ*(colD(:)-1), ...
                  'v',valD(:),'m',NZ*p,'n',NZ*q);
    return
end
Cout = sparse(rowMap(rowD), colMap(rowD) + NZ*(colD-1), valD, NZ*p, NZ*q);

end

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
%   layout   : (optional) 'separate' (default), 'packed' or 'triplets',     % MMP, 09/26/2026
%              selecting the output layout described below. May be given    % MMP, 09/26/2026
%              in place of Csize, since the two are told apart by type      % MMP, 09/26/2026
%
% OUTPUTS
%   C_gam_alp_beta : coefficient matrices, in one of three layouts          % MMP, 09/26/2026
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
% 'triplets' is 'separate' with each matrix M replaced by the struct        % MMP, 09/26/2026
%   struct('i',I,'j',J,'v',V,'m',size(M,1),'n',size(M,2))                   % MMP, 09/26/2026
% of column vectors, M == sparse(I,J,V,m,n), no repeated (I,J), no zero V,  % MMP, 09/26/2026
% in no particular order. The matrix is never formed: its g2*NR columns     % MMP, 09/26/2026
% carry the decision variables in 'copquadvar', so its column pointers      % MMP, 09/26/2026
% alone cost O(q) per cell to write and again to scan. Consumed by          % MMP, 09/26/2026
% 'lpis_sopvar/private/unpack_sheets'.                                      % MMP, 09/26/2026
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
% MMP, 09/30/2026: Split the 600-line primary (McCabe 69) into local
%                  functions along Sec. 6.1: semisep_setup (bases,
%                  reindexing, per-direction factor table),
%                  semisep_enumerate (the (beta,alpha,gamma) blocks) and
%                  semisep_assemble (memoized products, layouts); the
%                  primary keeps input handling. Code moved unchanged. The
%                  setup is a function of (G.Z, lims, g1) only and is now
%                  memoized (at most 16 setups, 16 MB), matched bit for
%                  bit on that key. Outputs bit-identical to the 09/29
%                  version (see the BEGIN/END block for the measurements).
% MMP, 09/29/2026: The Sec. 6.1 key table is now also written as the
%                  literal 11 x 3 matrix keyGBA, columns [gamma beta
%                  alpha], replacing str2double(num2cell(num2str(cllA')))
%                  (the column meaning was stated nowhere). The decimal
%                  cllA stays for the switch that keys on it. Removed the
%                  unused gamIdx (dec2base table, never read). Outputs
%                  bit-identical on all 3589 calls of the 1-D and 2-D
%                  container stability + Hinf builds and the 3-D heatNd
%                  build. Unprofiled against the HEAD copy on captured
%                  arguments: 3-D heatNd 2322 -> 2069 us per call
%                  (median), 9.20 -> 8.37 s per build; 1-D 336 -> 138 us.
%                  No change in the q scaling.
% MMP, 09/26/2026: New layout 'triplets' ('copquadvar' only): each output
%                  matrix is returned as its (i,j,v) triplets and never
%                  built. The NZ*p x NZ*q matrix of 'rearrangeCoef' has the
%                  decision variables on its q axis, so writing its column
%                  pointers (~108 MB per heavy 2-D Hinf call) and then
%                  scanning them in 'unpack_sheets' were 4.6 + 3.6 s of a
%                  19.7 s build, for ~4e3 nonzeros per call. Exact: with
%                  one beta and one alpha row each cell is ONE memoized
%                  'rearrangeCoef' result or empty, whose scatter is a
%                  bijection, so the triplets are that matrix's nonzeros;
%                  see 'unpack_sheets' for why its output depends only on
%                  their set. 'separate' and 'packed' are unchanged
%                  (operator products, 'sopquadvar', 'possopvar').
%                  Old -> new, one process, alternating, medians of 3,
%                  whole cx_exec builds with copquadvar's Zd change of the
%                  same date, programs isequal in every field: 2-D Hinf
%                  20.0 -> 10.8 s, its dual 20.1 -> 10.2 s, 2-D stability
%                  6.1 -> 4.1 s, and at degrees +1 (q = 5.5e5) 28.5 ->
%                  11.5 s; nv=3 poscopvar 2.8 -> 2.5 s, single space
%                  (q = 1.0e6) 46.9 -> 33.9 s; 1-D and the operator
%                  products unchanged. Peak memory in 'copquadvar' 216 ->
%                  17 MB (2-D Hinf), 1373 -> 13 MB (q = 5.5e5); cost per
%                  cell now O(nnz), no longer O(q*NZ). One exception: at
%                  nv=3 single space the cells are not wide and peak rises
%                  54 -> 69 MB (24 B/nonzero vs 16 B/nonzero + 8 B/column).
% MMP, 09/26/2026: 'rearrangeCoef' builds only the nonempty columns of
%                  X = reshape(C,NG,[]) when X is wide (p*q > 8*nnz(C)+1024).
%                  'copquadvar' puts the decision variables on q, so X had
%                  7.0e6 columns for 1.1e4 nonzeros, and forming X, the
%                  product and its reshape, each O(p*q), was 7.8 s of the
%                  routine's 12.8 s in the 2-D Hinf build. Now O(nnz(C)).
%                  Bit-identical: isequal on the 1-D and 2-D container
%                  programs and on every captured call. Old -> new, one
%                  process, medians: rearrangeCoef 14.0 -> 5.7 s (2-D Hinf),
%                  3.0 -> 1.3 s (2-D stability), 18.2 -> 7.8 s (2-D stability
%                  at degrees +1, q = 1.2e6); whole builds 26.3 -> 20.4 s,
%                  6.1 -> 4.2 s, 31.6 -> 21.2 s. 1-D builds and the operator
%                  products are unchanged: X is never wide there, so they keep
%                  the original lines. Transient memory per wide key drops from
%                  16*p*q bytes to O(nnz(C)) (108 -> 1.1 MB at the 2-D Hinf
%                  call). Peak memory and the remaining cost are both set by
%                  the output's NZ*q columns, which the interface fixes.
% MMP, 09/10/2026: Build the 'packed' outputs from accumulated triplets, one
%                  'sparse' call per gamma, instead of preallocating all-zero
%                  sparse matrices and filling them by subscripted
%                  assignment. Profiling a composition put 67.5% of this
%                  routine in that single assignment statement at two
%                  pass-through variables and degree 3, and 79.5% at three;
%                  its cost per nonzero moved rose from 67.9 to 131.0 ns as
%                  the problem grew, which is the signature of rebuilding
%                  the sparse structure on every write. CLAUDE.md's rule --
%                  build sparse matrices once from triplets, never
%                  index-assign in a loop -- applied literally. Only the
%                  'packed' layout is affected; the 'separate' layout
%                  already assigned whole cells.
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
        || ~any(strcmpi(layout,{'separate','packed','triplets'}))           % MMP, 09/26/2026
%       || ~any(strcmpi(layout,{'separate','packed'}))                       % MMP, 08/30/2026 % MMP, 09/26/2026 (was)
    error(['int_semisep: layout must be ''separate'', ''packed'' or ' ...
           '''triplets''.']);                                               % MMP, 09/26/2026
%   error('int_semisep: layout must be ''separate'' or ''packed''.');        % MMP, 08/30/2026 % MMP, 09/26/2026 (was)
end                                                                          % MMP, 08/30/2026
packed = strcmpi(layout,'packed');                                           % MMP, 08/30/2026
trip   = strcmpi(layout,'triplets');    % 'separate' shape, triplet entries % MMP, 09/26/2026

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
        if trip     % CG itself is the output; hand over its nonzeros       % MMP, 09/26/2026
            [ti,tj,tv] = find(CG);                                          % MMP, 09/26/2026
            CG = struct('i',ti(:),'j',tj(:),'v',tv(:), ...
                        'm',size(CG,1),'n',size(CG,2));                     % MMP, 09/26/2026
        end                                                                 % MMP, 09/26/2026
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

% BEGIN MMP, 09/30/2026: split into setup / enumerate / assemble, setup
% memoized. What was the rest of this function is now three local
% functions (below), one per object of Sec. 6.1 of
% sopvar_implementation_notes.pdf:
%   semisep_setup     the bases ZL, ZR, the reindexing of 'rearrangeCoef'
%                     and the per-direction factors of the eleven keys;
%   semisep_enumerate for each requested (beta,alpha) row pair, the gammas
%                     with p(gamma,alpha,beta) nonzero and their factors;
%   semisep_assemble  each block's coefficients (Sec. 6.2) in the layout.
% Lines moved unchanged; only the lines marked 09/30 are new. Nothing was
% deleted apart from blank lines. The setup depends on (G.Z, lims, g1)
% alone and the callers repeat it, so it is fetched from a bounded memo
% (memo_setup, which states the measurements).
S = memo_setup(ZG,lims,g1);                                                 % MMP, 09/30/2026
ZL = S.ZL;      ZR = S.ZR;                                                  % MMP, 09/30/2026
[blk,kseqs] = semisep_enumerate(S,idxbeta,idxalpha);                        % MMP, 09/30/2026
C_gam_alp_beta = semisep_assemble(S,blk,kseqs,CG,NG,g1,g2,nbeta,nalpha, ...
                                  packed,trip);                             % MMP, 09/30/2026
% END MMP, 09/30/2026

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = setup_key(ZG,lims,g1)                                        % MMP, 09/30/2026
% The setup key (ZG, lims, g1) as the 64-bit patterns of g1, the counts and % MMP, 09/30/2026
% every entry, so equal keys give bit-identical setups; isequal on the      % MMP, 09/30/2026
% values would equate 0 and -0. [] (never cached) unless lims and every     % MMP, 09/30/2026
% ZG{i} are real full double and each ZG{i} a column: those fix the class   % MMP, 09/30/2026
% and shape of the setup, and the counts then fix every size. lims is       % MMP, 09/30/2026
% numel(ZG) x 2 (checked above). A few builtins, independent of numel(ZG):  % MMP, 09/30/2026
% an entry-by-entry comparison cost 6-22 us per call (measured).            % MMP, 09/30/2026
key = [];                                                                   % MMP, 09/30/2026
if isa(lims,'double') && ~issparse(lims) && isreal(lims) ...
        && all(cellfun('isclass',ZG,'double')) && all(cellfun('isreal',ZG)) ...
        && all(cellfun('size',ZG,2)==1)                                     % MMP, 09/30/2026
    v = vertcat(ZG{:});         % sparse if any ZG{i} is                    % MMP, 09/30/2026
    if ~issparse(v)                                                         % MMP, 09/30/2026
        key = typecast([g1; numel(ZG); cellfun('prodofsize',ZG(:)); ...
                        lims(:); v],'uint64');                              % MMP, 09/30/2026
    end                                                                     % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026


function S = memo_setup(ZG,lims,g1)                                         % MMP, 09/30/2026
% The setup of (ZG, lims, g1), from the memo when it holds that key.        % MMP, 09/30/2026
% copquadvar calls int_semisep once per (i,j) basis pair with G.Z and the   % MMP, 09/30/2026
% domain fixed: 3251 calls in the 3-D heatNd build, 9 distinct keys, 61     % MMP, 09/30/2026
% calls whose key differs from the previous call's, 3.5 s of setup. The     % MMP, 09/30/2026
% 1-D and 2-D builds have 9 and 10 keys but 68 and 54 changes, so the memo  % MMP, 09/30/2026
% holds several setups, not the last one only. Keys match bit for bit       % MMP, 09/30/2026
% (setup_key); the last match is tried first. Bounded: at most 16 setups    % MMP, 09/30/2026
% and 2^20 rowMap entries in all (rowMap and colMap have g1*NL*NR entries   % MMP, 09/30/2026
% each: 16 MB), oldest dropped first. A setup is of spatial size and holds  % MMP, 09/30/2026
% no G.C, so nothing q-sized is kept.                                       % MMP, 09/30/2026
persistent keys Ss last                                                     % MMP, 09/30/2026
key = setup_key(ZG,lims,g1);                                                % MMP, 09/30/2026
if ~isempty(key) && ~isempty(last) && isequal(key,keys{last})               % MMP, 09/30/2026
    S = Ss{last};                                                           % MMP, 09/30/2026
    return                                                                  % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026
if ~isempty(key)                                                            % MMP, 09/30/2026
    for i = 1:numel(keys)                                                   % MMP, 09/30/2026
        if isequal(key,keys{i}),    S = Ss{i};  last = i;   return,     end % MMP, 09/30/2026
    end                                                                     % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026
S = semisep_setup(ZG,lims,g1);                                              % MMP, 09/30/2026
if ~isempty(key) && numel(S.rowMap)<=2^20                                   % MMP, 09/30/2026
    keys{end+1} = key;      Ss{end+1} = S;                                  % MMP, 09/30/2026
    nmap = cellfun(@(s) numel(s.rowMap),Ss);                                % MMP, 09/30/2026
    while numel(Ss)>16 || sum(nmap)>2^20                                    % MMP, 09/30/2026
        keys(1) = [];   Ss(1) = [];     nmap(1) = [];                       % MMP, 09/30/2026
    end                                                                     % MMP, 09/30/2026
    last = numel(Ss);                                                       % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = semisep_setup(ZG,lims,g1)                                      % MMP, 09/30/2026
% Everything of int_semisep that depends on the key (ZG, lims, g1) alone,   % MMP, 09/30/2026
% Sec. 6.1: the common bases ZL = ZR, the permutation of 'rearrangeCoef'    % MMP, 09/30/2026
% (rowMap, colMap), the factor Ci_key_ell{key,ell} of each key and          % MMP, 09/30/2026
% direction with the classes of equal factors (clsOf, nCls, pwCls), and     % MMP, 09/30/2026
% the key lists key_of{beta,alpha}. Nothing here reads G.C, the index rows  % MMP, 09/30/2026
% or the layout, which is what makes caching S exact.                       % MMP, 09/30/2026
ns3a = numel(ZG);                                                           % MMP, 09/30/2026
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
% gamIdx = fliplr(dec2base(0:3^ns3a-1,3,ns3a)-'0') + 1;                     % MMP, 09/29/2026 (was)
% Removed: never read; gam_lin in the key loop encodes gamma in base 3.     % MMP, 09/29/2026

% The eleven (gamma,beta,alpha) keys for which p(gamma,alpha,beta) is
% nonzero, i.e. the eleven cases of the integral table in Sec. 6.1.
cllA = [111   212   313   221   331   222   232   332   223   323   333];
% The same eleven keys as rows [gamma beta alpha], row k the digits of      % MMP, 09/29/2026
% cllA(k): column 1 gamma, 2 beta, 3 alpha. Replaces cllA after the switch. % MMP, 09/29/2026
keyGBA = [1 1 1; 2 1 2; 3 1 3; 2 2 1; 3 3 1; 2 2 2; ...
          2 3 2; 3 3 2; 2 2 3; 3 2 3; 3 3 3];                               % MMP, 09/29/2026

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

% cllA = str2double(num2cell(num2str(cllA')));                              % MMP, 09/29/2026 (was)
% The literal above, not a number->string->number parse (165 us per call).  % MMP, 09/29/2026
cllA = keyGBA;                                                              % MMP, 09/29/2026

% The eleven keys share only EIGHT distinct per-direction factors: the switch  % MMP, 09/10/2026
% above handles {111,212,313} in one branch and {221,331} in another, so       % MMP, 09/10/2026
% Ci_key_ell takes 8 values per direction and the Kronecker product Csep       % MMP, 09/10/2026
% below takes 8^ns3a, not 11^ns3a. Every repeat costs a redundant              % MMP, 09/10/2026
% 'rearrangeCoef', which is 11^ns3a/8^ns3a = 1.89x at two pass-through         % MMP, 09/10/2026
% variables, 2.60x at three and 3.57x at four. The duplicates are found by     % MMP, 09/10/2026
% COMPARING the built factors rather than by hard-coding the branch labels,    % MMP, 09/10/2026
% so the grouping stays correct if the table or the domains change; the        % MMP, 09/10/2026
% comparison is 11^2 per direction on matrices of a few dozen nonzeros.        % MMP, 09/10/2026
nkey  = size(cllA,1);                                                        % MMP, 09/10/2026
repOf = repmat((1:nkey).',1,ns3a);                                           % MMP, 09/10/2026
for ell = 1:ns3a                                                             % MMP, 09/10/2026
    for k = 2:nkey                                                           % MMP, 09/10/2026
        for j = 1:k-1                                                        % MMP, 09/10/2026
            if isequal(Ci_key_ell{k,ell},Ci_key_ell{j,ell})                   % MMP, 09/10/2026
                repOf(k,ell) = repOf(j,ell);                                 % MMP, 09/10/2026
                break                                                        % MMP, 09/10/2026
            end                                                              % MMP, 09/10/2026
        end                                                                  % MMP, 09/10/2026
    end                                                                      % MMP, 09/10/2026
end                                                                          % MMP, 09/10/2026
% Renumber the representatives densely per direction, so the memo is indexed  % MMP, 09/10/2026
% in mixed radix over the DISTINCT counts (8 per direction in the current      % MMP, 09/10/2026
% table) rather than over 11, and holds prod(nCls) slots instead of           % MMP, 09/10/2026
% 11^ns3a.                                                                    % MMP, 09/10/2026
clsOf = zeros(nkey,ns3a);   nCls = zeros(1,ns3a);                            % MMP, 09/10/2026
for ell = 1:ns3a                                                             % MMP, 09/10/2026
    [~,~,ic] = unique(repOf(:,ell));                                         % MMP, 09/10/2026
    clsOf(:,ell) = ic;                                                       % MMP, 09/10/2026
    nCls(ell)    = max(ic);                                                  % MMP, 09/10/2026
end                                                                          % MMP, 09/10/2026
pwCls = cumprod([1,nCls(1:end-1)]);                                          % MMP, 09/10/2026

pw3 = 3.^(0:ns3a-1)';                                                        % MMP, 08/30/2026
key_of = cell(3,3);                                                          % MMP, 08/30/2026
for k = 1:size(cllA,1)                                                       % MMP, 08/30/2026
    key_of{cllA(k,2),cllA(k,3)} = [key_of{cllA(k,2),cllA(k,3)}, k];          % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026
S = struct('ZL',{ZL},'ZR',{ZR},'NL',NL,'NR',NR,'rowMap',rowMap, ...
           'colMap',colMap,'Ci_key_ell',{Ci_key_ell},'cllA',cllA, ...
           'clsOf',clsOf,'nCls',nCls,'pwCls',pwCls,'pw3',pw3, ...
           'key_of',{key_of});                                              % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [blk,kseqs] = semisep_enumerate(S,idxbeta,idxalpha)                % MMP, 09/30/2026
% The blocks to build, Sec. 6.1: p(gamma,alpha,beta) is the product over    % MMP, 09/30/2026
% directions of the per-direction table, so for a (beta,alpha) row pair     % MMP, 09/30/2026
% the gammas are the Cartesian product of the per-direction key options.    % MMP, 09/30/2026
% Row t of blk is [gamma cell, beta row, alpha row, memo slot] of block t,  % MMP, 09/30/2026
% in the order the loop below met it; kseqs(t,:) its key per direction.     % MMP, 09/30/2026
% Index work only, on at most nbeta*nalpha*2^ns3a rows.                     % MMP, 09/30/2026
ns3a = numel(S.ZL);                                                         % MMP, 09/30/2026
nbeta = size(idxbeta,1);        nalpha = size(idxalpha,1);                  % MMP, 09/30/2026
key_of = S.key_of;  cllA = S.cllA;  clsOf = S.clsOf;                        % MMP, 09/30/2026
pwCls = S.pwCls;    pw3 = S.pw3;                                            % MMP, 09/30/2026
% A (beta,alpha) pair admits at most max numel(key_of) keys per direction.  % MMP, 09/30/2026
nmax = nbeta*nalpha*max(cellfun('prodofsize',key_of(:)))^ns3a;              % MMP, 09/30/2026
blk = zeros(nmax,4);    kseqs = zeros(nmax,ns3a);   nb = 0;                 % MMP, 09/30/2026

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
    % Walk the mixed radix for the keys, gamma, and the memo slot. Csep is    % MMP, 09/10/2026
    % NOT built here: it is only needed on a memo miss, so the Kronecker      % MMP, 09/10/2026
    % chain is skipped on a hit along with the rearrangeCoef.                 % MMP, 09/10/2026
    gam_lin = 1;    rem_c = comb;    memoLin = 1;                            % MMP, 09/10/2026
    kseq = zeros(1,ns3a);                                                    % MMP, 09/10/2026
    for ell = 1:ns3a                                                         % MMP, 08/30/2026
        sel = mod(rem_c,nopt(ell)) + 1;                                      % MMP, 08/30/2026
        rem_c = floor(rem_c/nopt(ell));                                      % MMP, 08/30/2026
        key_idx = opt{ell}(sel);                                             % MMP, 08/30/2026
        kseq(ell) = key_idx;                                                 % MMP, 09/10/2026
        memoLin = memoLin + (clsOf(key_idx,ell)-1)*pwCls(ell);               % MMP, 09/10/2026
%       Csep = kron(Csep,Ci_key_ell{key_idx,ell});                           % MMP, 09/10/2026 (was)
        gam_lin = gam_lin + (cllA(key_idx,1)-1)*pw3(ell);                    % MMP, 08/30/2026
    end                                                                      % MMP, 08/30/2026
    un_indices_gamma = gam_lin;                                              % MMP, 08/30/2026
    % Record the block; semisep_assemble builds it.                         % MMP, 09/30/2026
    nb = nb + 1;                                                            % MMP, 09/30/2026
    blk(nb,:) = [un_indices_gamma,un_indices_beta, ...
                 un_indices_alpha,memoLin];                                 % MMP, 09/30/2026
    kseqs(nb,:) = kseq;                                                     % MMP, 09/30/2026
end                                                                          % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026
end                                                                          % MMP, 08/30/2026

% Repeated rows in IDXBETA or IDXALPHA need no special handling: the loop     % MMP, 08/30/2026
% above is driven by the row indices themselves, so every occurrence is       % MMP, 08/30/2026
% filled. The earlier version resolved a multi-index back to its first        % MMP, 08/30/2026
% matching row and needed a separate pass to copy the duplicates.             % MMP, 08/30/2026
blk = blk(1:nb,:);      kseqs = kseqs(1:nb,:);                              % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C_gam_alp_beta = semisep_assemble(S,blk,kseqs,CG,NG,g1,g2, ...
                                           nbeta,nalpha,packed,trip)        % MMP, 09/30/2026
% The coefficients of every block of blk (Sec. 6.2): the Kronecker product  % MMP, 09/30/2026
% of its per-direction factors kseqs, applied to CG and re-indexed onto ZL, % MMP, 09/30/2026
% ZR by 'rearrangeCoef', memoized per slot since blocks share factors;      % MMP, 09/30/2026
% written into the requested layout. Blocks are taken in blk's order, so    % MMP, 09/30/2026
% each slot is built from the same block, and 'rearrangeCoef' runs in the   % MMP, 09/30/2026
% same order, as when this was one loop.                                    % MMP, 09/30/2026
ns3a = numel(S.ZL);     NL = S.NL;      NR = S.NR;      nCls = S.nCls;      % MMP, 09/30/2026
Ci_key_ell = S.Ci_key_ell;  rowMap = S.rowMap;  colMap = S.colMap;          % MMP, 09/30/2026
% 'packed' needs 3^ns3a matrices rather than 3^ns3a*nbeta*nalpha of them.     % MMP, 08/30/2026
if packed                                                                    % MMP, 08/30/2026
    % The packed outputs are accumulated as triplets and built with one      % MMP, 09/10/2026
    % 'sparse' call each, below the loop. They used to be preallocated as    % MMP, 09/10/2026
    % all-zero sparse matrices and filled by subscripted assignment, which   % MMP, 09/10/2026
    % rebuilds the whole sparse structure on every one of the 11^ns3a        % MMP, 09/10/2026
    % writes: that single statement measured 67.5% of this routine at two    % MMP, 09/10/2026
    % pass-through variables and degree 3, rising to 79.5% at three, and     % MMP, 09/10/2026
    % its cost per nonzero moved GREW with the problem (67.9 -> 131.0 ns).   % MMP, 09/10/2026
    % One slot per (gamma,beta,alpha): the eleven rows of 'cllA' are         % MMP, 09/10/2026
    % distinct (gamma,beta,alpha) triples, so for a fixed (beta,alpha) the   % MMP, 09/10/2026
    % gammas differ, and gamma-tuple -> gam_lin is a bijection. Each block   % MMP, 09/10/2026
    % is therefore written at most once and nothing needs summing.          % MMP, 09/10/2026
    nGam = 3^ns3a;                                                           % MMP, 09/10/2026
    accI = cell(nGam,nbeta*nalpha);                                          % MMP, 09/10/2026
    accJ = cell(nGam,nbeta*nalpha);                                          % MMP, 09/10/2026
    accV = cell(nGam,nbeta*nalpha);                                          % MMP, 09/10/2026
%   C_gam_alp_beta = cell(3^ns3a,1);                                         % MMP, 09/10/2026 (was)
%   for idx_C = 1:numel(C_gam_alp_beta)                                      % MMP, 09/10/2026 (was)
%       C_gam_alp_beta{idx_C} = sparse(g1*NL*nbeta,g2*NR*nalpha);            % MMP, 09/10/2026 (was)
%   end                                                                      % MMP, 09/10/2026 (was)
elseif trip                                                                 % MMP, 09/26/2026
    % Unwritten cells: no triplets, the dimensions of the all-zero matrix.  % MMP, 09/26/2026
    C_gam_alp_beta = cell(3^ns3a,nbeta,nalpha);                             % MMP, 09/26/2026
    C_gam_alp_beta(:) = {struct('i',zeros(0,1),'j',zeros(0,1), ...
        'v',zeros(0,1),'m',g1*NL,'n',g2*NR)};                               % MMP, 09/26/2026
else                                                                         % MMP, 08/30/2026
    C_gam_alp_beta = cell(3^ns3a,nbeta,nalpha);
    for idx_C = 1:numel(C_gam_alp_beta)
        C_gam_alp_beta{idx_C} = sparse(g1*NL,g2*NR);
    end
end                                                                          % MMP, 08/30/2026
memoC = cell(prod([nCls,1]),1);                                              % MMP, 09/10/2026
for t = 1:size(blk,1)                                                       % MMP, 09/30/2026
    un_indices_gamma = blk(t,1);        un_indices_beta = blk(t,2);         % MMP, 09/30/2026
    un_indices_alpha = blk(t,3);        memoLin = blk(t,4);                 % MMP, 09/30/2026
    kseq = kseqs(t,:);                                                      % MMP, 09/30/2026

    % Convert
    %
    %   (I_g1 \otimes Zmix')*(I_g1 \otimes Csep')*CG
    %
    % into
    %
    %   (I_g1 \otimes ZL')*Cnew*(I_g2 \otimes ZR).
    % Keyed on the representative tuple, so the 11^ns3a enumerated blocks     % MMP, 09/10/2026
    % share the 8^ns3a distinct results. The size checks stay on the miss     % MMP, 09/10/2026
    % path, so every distinct Csep and Cnew is still validated once.          % MMP, 09/10/2026
    if isempty(memoC{memoLin})                                               % MMP, 09/10/2026
        Csep = 1;                                                            % MMP, 09/10/2026
        for ell = 1:ns3a                                                     % MMP, 09/10/2026
            Csep = kron(Csep,Ci_key_ell{kseq(ell),ell});                     % MMP, 09/10/2026
        end                                                                  % MMP, 09/10/2026
        if ~isequal(size(Csep),[NG,NL*NR])
            error('int_semisep: internal Csep size mismatch.');
        end
%       Cmiss = rearrangeCoef(Csep,CG,g1,g2,rowMap,colMap,NL);              % MMP, 09/10/2026 % MMP, 09/26/2026 (was)
%       if ~isequal(size(Cmiss),[g1*NL,g2*NR])                               % MMP, 09/10/2026 % MMP, 09/26/2026 (was)
        % 'triplets': Cmiss is a triplet struct and its size is its m, n.   % MMP, 09/26/2026
        Cmiss = rearrangeCoef(Csep,CG,g1,g2,rowMap,colMap,NL,trip);         % MMP, 09/26/2026
        if trip,  szC = [Cmiss.m,Cmiss.n];  else,  szC = size(Cmiss);  end  % MMP, 09/26/2026
        if ~isequal(szC,[g1*NL,g2*NR])                                      % MMP, 09/26/2026
            error('int_semisep: output coefficient size mismatch.');
        end
        memoC{memoLin} = Cmiss;                                              % MMP, 09/10/2026
    end                                                                      % MMP, 09/10/2026
    Cnew = memoC{memoLin};                                                   % MMP, 09/10/2026

    if packed                                                                % MMP, 08/30/2026
        % Stash this block's nonzeros, shifted onto its (beta,alpha) offset. % MMP, 09/10/2026
        [ri,ci,vi] = find(Cnew);                                             % MMP, 09/10/2026
        slot = un_indices_beta + nbeta*(un_indices_alpha - 1);               % MMP, 09/10/2026
        accI{un_indices_gamma,slot} = ri + g1*NL*(un_indices_beta - 1);      % MMP, 09/10/2026
        accJ{un_indices_gamma,slot} = ci + g2*NR*(un_indices_alpha - 1);     % MMP, 09/10/2026
        accV{un_indices_gamma,slot} = vi;                                    % MMP, 09/10/2026
%       rows = (1:(g1*NL)) + g1*NL*(un_indices_beta - 1);                    % MMP, 09/10/2026 (was)
%       cols = (1:(g2*NR)) + g2*NR*(un_indices_alpha - 1);                   % MMP, 09/10/2026 (was)
%       C_gam_alp_beta{un_indices_gamma}(rows,cols) = Cnew;                  % MMP, 09/10/2026 (was)
    else                                                                     % MMP, 08/30/2026
        ind_Cgam_alpha = sub2ind(size(C_gam_alp_beta), un_indices_gamma, un_indices_beta, un_indices_alpha);
        C_gam_alp_beta{ind_Cgam_alpha} = Cnew;
    end                                                                      % MMP, 08/30/2026
end                                                                         % MMP, 09/30/2026

% One 'sparse' per gamma from the accumulated triplets. A gamma that was      % MMP, 09/10/2026
% never written concatenates to empty, which 'sparse' turns into the          % MMP, 09/10/2026
% all-zero matrix the preallocation used to supply.                          % MMP, 09/10/2026
if packed                                                                    % MMP, 09/10/2026
    C_gam_alp_beta = cell(nGam,1);                                           % MMP, 09/10/2026
    for k = 1:nGam                                                           % MMP, 09/10/2026
        C_gam_alp_beta{k} = sparse(vertcat(accI{k,:}),vertcat(accJ{k,:}), ...% MMP, 09/10/2026
            vertcat(accV{k,:}),g1*NL*nbeta,g2*NR*nalpha);                    % MMP, 09/10/2026
    end                                                                      % MMP, 09/10/2026
end                                                                          % MMP, 09/10/2026
end                                                                         % MMP, 09/30/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Cout = rearrangeCoef(Csep,C,p,q,rowMap,colMap,NZ,trip)             % MMP, 09/26/2026
% function Cout = rearrangeCoef(Csep,C,p,q,rowMap,colMap,NZ)                 % MMP, 08/30/2026 % MMP, 09/26/2026 (was)
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
%   When that reshape is wide -- q carries the decision variables and its   % MMP, 09/26/2026
%   p*q columns are nearly all empty -- only its nonempty columns are       % MMP, 09/26/2026
%   formed, from C's triplets.                                              % MMP, 09/26/2026
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
%   trip   : (optional, default false) return Cout as the triplet struct    % MMP, 09/26/2026
%            of the 'triplets' layout instead of building it                % MMP, 09/26/2026
%
% OUTPUT
%   Cout   : NZ*p by NZ*q sparse matrix
if nargin<8,    trip = false;   end                                         % MMP, 09/26/2026

% BEGIN MMP, 09/26/2026: wide X. The line after END forms X = reshape(C,NG,[])
% with p*q columns, the product D0 = Csep.'*X over all of them, and a reshape
% of D0: each O(p*q), whatever nnz(C). In 'copquadvar' q carries the
% decision variables and nearly every column of X is empty (2-D Hinf:
% p*q = 7.0e6 against nnz(C) = 1.1e4), which made that line 7.8 s of this
% routine's 12.8 s. Here the nonempty columns of X alone are built from C's
% triplets, O(nnz(C)) plus one read of C, and mapped back.
%   Exact: sparse*sparse computes each column of D0 from the same column of
% X alone, so every entry is the same sum of the same products in the same
% order; the result is bit-identical (isequal on the 1-D and 2-D container
% Hinf and stability programs, and on every call captured from them and
% from the operator products).
%   The rule is a measured crossover over synthetic shapes at p = 8, 25, 64:
% under it no shape ran slower than the original. Narrow X, e.g. every call
% from the operator products, where p*q is about nnz(C), keeps the original.
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
    if trip     % the scatter's triplets: bijective, so no duplicates       % MMP, 09/26/2026
        Cout = struct('i',rowMap(rowD(:)),'j',colMap(rowD(:))+NZ*(colD(:)-1), ...
                      'v',valD(:),'m',NZ*p,'n',NZ*q);                       % MMP, 09/26/2026
        return                                                              % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
    Cout = sparse(rowMap(rowD), colMap(rowD) + NZ*(colD-1), valD, NZ*p, NZ*q);
    return
end
% END MMP, 09/26/2026

D = reshape(sparse(Csep).' * reshape(sparse(C),size(Csep,1),[]), [], q);

if isempty(rowMap)
    Cout = D;
    return
end
if size(D,1)~=numel(rowMap)
    error('rearrangeCoef: input dimensions are inconsistent.');
end

[rowD,colD,valD] = find(D);
if trip         % as in the wide branch above                               % MMP, 09/26/2026
    Cout = struct('i',rowMap(rowD(:)),'j',colMap(rowD(:))+NZ*(colD(:)-1), ...
                  'v',valD(:),'m',NZ*p,'n',NZ*q);                           % MMP, 09/26/2026
    return                                                                  % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
Cout = sparse(rowMap(rowD), colMap(rowD) + NZ*(colD-1), valD, NZ*p, NZ*q);

end

function [prog,M,info] = posmult_cdopvar(prog,dims,spaces,dom,w,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,M,INFO] = POSMULT_CDOPVAR(PROG,DIMS,SPACES,DOM,W,OPTIONS) declares
% the positive semidefinite, self-adjoint 'cdopvar' decision operator
%
%       M = E* Phi E,
%       Phi(theta) = sum_c g_c(theta) L_c(theta)' Q_c L_c(theta),
%
% on X = X_1 x ... x X_Ms, X_k = L_2^{m_k}[s^k] (R^{m_k} when s^k is
% empty), with one Gram matrix Q_c >= 0 per entry c of OPTIONS.psatz. Here
% theta = (theta_1,...,theta_nv) ranges over the box of the registry, the
% sorted union of the s^k; (E x)(theta) = [x_1(theta^1); ...;
% x_Ms(theta^Ms)], theta^k the entries of theta of the variables of s^k
% (x_k itself on R^{m_k}), and E* integrates the other entries over the
% box; L_c(theta) = blkdiag_k(I_{m_k} kron Z_{c,k}(theta)), Z_{c,k} the
% vector of the monomials theta^p with p_d <= max(w_k(d) - o_c, 0) for
% every d, w_k the weight degrees of space k (W) and o_c the degree
% reduction of term c (OPTIONS.psatz_offset); and g_c, the weight of code
% c, is 1 (code 0), prod_d (theta_d-a_d)(b_d-theta_d) (code 1),
% (theta_d-a_d)/(b_d-a_d) (code 2d+1) or (b_d-theta_d)/(b_d-a_d) (code
% 2d+2). Every g_c is nonnegative on the box, so Phi(theta) >= 0 there and
% M >= 0.
%
% INPUT
% - prog:    'struct', the LPI program to which the Grams Q_c and the
%            constraints Q_c >= 0 are added;
% - dims:    'double', Ms x 1 array of positive integers, dims(k) = m_k, or
%            a scalar for every space;
% - spaces:  1 x Ms 'cell', spaces{k} the 'cellstr' of the variable names
%            of space k, empty for R^{m_k}; a 'cellstr' is one space;
% - dom:     'double', nv x 2 array, row d the interval [a_d,b_d] of
%            registry variable d; or a 1 x 2 row for every variable; or a
%            'struct' with fields vars (a 'cellstr' of the names) and dom
%            (one row per name);
% - w:       (optional) 'double', the weight degrees, nonnegative integers:
%            a scalar or a 1 x nv row over the registry, for every space;
%            or a 1 x Ms 'cell' with one such entry per space (default 1);
% - options: (optional) 'struct' with the optional fields
%              psatz         'double' row of codes, one term per entry,
%                            repeated codes allowed: 0, 1, 2d+1 or 2d+2,
%                            d = 1..nv over the registry (default 0);
%              psatz_offset  'double', the degree reduction o_c of each
%                            term: a scalar for every term or one entry per
%                            code (default 1 for code 1, 0 otherwise).
%
% OUTPUT
% - prog:    'struct', PROG with the Grams declared in the order of
%            OPTIONS.psatz, each with its constraint Q_c >= 0;
% - M:       'cdopvar' object, Ms x Ms, on X; M = M* >= 0;
% - info:    'struct' with fields
%              psatz   1 x nt 'double', the codes of the nt terms;
%              offset  1 x nt 'double', the degree reductions o_c;
%              deg     1 x nt 'cell', deg{c} the 1 x Ms 'cell' of 'struct'
%                      with field int, the caps of Z_{c,k};
%              Qcell   1 x nt 'cell', Qcell{c} the Ms x Ms 'cell' of
%                      'cellstr' arrays, the names of the entries of the
%                      blocks (k,l) of Q_c;
%              gram    1 x nt 'double', the dimension of each Q_c;
%              ndec    'double', the number of decision variables of M.
%
% See also POSCOPVAR_LIFT, LIFT_COPVAR, POSCOPVAR, COPQUADVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - posmult_cdopvar
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
% Initial coding MMP, 10/07/2026. The positive multiplier of the separated
%                form A'*M*A proposed from the proof program: the parameter
%                of a positive operator is the pointwise PSD polynomial
%                matrix M, parameterised as a sum of matrix SOS terms with
%                the interval Positivstellensatz weights, each term one
%                degree-counted Gram.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers of
%                the local functions; code unchanged. The deleted header
%                text was the Positivstellensatz background of the weights
%                with its notes references, the role of M in the proof
%                program, the comparison with copquadvar, and the cost and
%                measurement notes.

% Four inputs required; W = 1 and OPTIONS = struct() when omitted
if nargin<4
    error("Not enough input arguments.")
end
if nargin<5 || isempty(w),          w = 1;              end
if nargin<6 || isempty(options),    options = struct(); end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end

% % % Spaces, for the per-space degrees and the all-multiplier 'include'.
[meta,sp] = parse_copvar_spaces(dims,spaces,dom);
Ms = numel(sp);     nv = numel(meta.vars);
isR = cellfun(@isempty,sp);

% % % Weight degrees, one 1 x nv row per space.
if iscell(w)
    if numel(w)~=Ms
        error("A cell 'w' should have one entry per space.")
    end
    wsp = cellfun(@(v) expand_row(v,nv),w,'UniformOutput',false);
else
    wsp = repmat({expand_row(w,nv)},1,Ms);
end

% % % Terms
% The codes, a row (0 when absent), each 0, 1 or a face code
codes = 0;
if isfield(options,'psatz') && ~isempty(options.psatz)
    codes = reshape(double(options.psatz),1,[]);
end
if any(~ismember(codes,[0,1,3:2*nv+2]))
    error("'psatz' entries should be 0, 1, or face codes 2d+1 / 2d+2 with d in 1..nv.")
end
% Degree reduction per term: 1 for code 1 and 0 otherwise, or psatz_offset,
% a scalar or one entry per term
nt = numel(codes);
off = double(codes==1);                 % deg S1 = deg S0 - 2: one per direction
if isfield(options,'psatz_offset') && ~isempty(options.psatz_offset)
    off = reshape(double(options.psatz_offset),1,[]);
    if isscalar(off),   off = repmat(off,1,nt);
    elseif numel(off)~=nt
        error("'psatz_offset' should be a scalar or have one entry per psatz term.")
    end
end

% The multiplier component of every L2 space is row 1 of the multi-index
% enumeration (all ones, whatever 'sep'); an R space has only its identity.
incl = cell(1,Ms);
for k = 1:Ms
    if ~isR(k),     incl{k} = 1;    end
end

% One Gram per term at the reduced degrees, the terms summed into M
M = [];
info = struct('psatz',codes,'offset',off,'deg',{cell(1,nt)},...
              'Qcell',{cell(1,nt)},'gram',zeros(1,nt),'ndec',0);
for c = 1:nt
    % Caps of term c on every space: the weight degrees less the reduction
    deg = cell(1,Ms);
    for k = 1:Ms
        deg{k} = struct('int',max(wsp{k}-off(c),0));
    end
    % Term c: its Gram on the multiplier basis operators, weight of code c
    popt = struct('psatz',codes(c),'include',{incl});
    [prog,Mc,Qc] = poscopvar(prog,dims,spaces,dom,deg,popt);
    % Record the term; its Gram dimension; add it to M
    info.deg{c} = deg;      info.Qcell{c} = Qc;
    info.gram(c) = sum(arrayfun(@(i) size(Qc{i,i},1),1:size(Qc,1)));
    if isempty(M),  M = Mc;     else,   M = M + Mc;  end
end

% % % Prune every block to its occupied monomials. 'copquadvar' returns a
% block on the bases its pair integrals can reach, and a multiplier-only
% block occupies right degree 0 and left degree 2w of them only (measured
% 10/07/2026: 3-D, w = 1, 8000 x 8000 coefficient grids per cell with
% 421875 of 64e6 columns in use). The composition A'*M*A by
% '@cdopvar/mtimes' carries M's bases through every intermediate, so the
% padding multiplies its time and memory. The operator is unchanged.
[Mo,Mi] = size(M.C);
for i = 1:Mo
    for j = 1:Mi
        if isa(M.C{i,j},'sdopvar'),     M.C{i,j} = prune_block(M.C{i,j});   end
    end
end
% Number of decision variables of M
info.ndec = numel(M.Zd);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = expand_row(v,nv)
% V = EXPAND_ROW(V,NV) returns the weight degrees V as a row with one entry
% per registry variable.
% INPUT  v:  'double', nonnegative integers: a scalar, or nv entries;
%        nv: 'double', the number of registry variables.
% OUTPUT v:  1 x nv 'double' (1 x 1 for a scalar V when nv = 0).

% Nonnegative integers; a scalar repeated, else one entry per variable
if ~isnumeric(v) || any(v(:)<0) || any(v(:)~=round(v(:)))
    error("Weight degrees should be nonnegative integers.")
end
v = reshape(v,1,[]);
if isscalar(v)
    v = repmat(v,1,max(nv,1));
elseif numel(v)~=nv
    error("A weight degree should be a scalar or have one entry per registry variable.")
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = prune_block(B)
% B = PRUNE_BLOCK(B) removes from the monomial bases of the block B, per
% direction, the degrees that carry no coefficient: a degree is dropped
% only when neither A nor any row of B is nonzero there, so the operator
% is unchanged, and the canonical multiplier form survives (right degree
% 0 stays occupied wherever a multiplier cell has content).
% INPUT  B: 'sdopvar' object.
% OUTPUT B: 'sdopvar' object, the same operator on the pruned bases (B
%           itself when no degree is empty).

% Coefficient sizes: rows (matrix row, left monomial), columns (matrix
% column, right monomial), monomial fastest
m = B.dims(1);  n = B.dims(2);
ZL = B.ZL;      ZR = B.ZR;
NL = prod([cellfun(@numel,ZL),1]);  NR = prod([cellfun(@numel,ZR),1]);
nR = m*NL;      nC = n*NR;
pA = B.params.A;    pB = B.params.B;
% Occupied left and right monomials: those of the nonzeros of A and of the
% columns of B with a nonzero, in every cell stored at full size
occL = false(NL,1); occR = false(NR,1);
for ii = 1:numel(pA)
    if numel(pA{ii})~=nR*nC,    continue,   end
    lin = find(pA{ii});
    if ~isempty(pB{ii}),    [~,cB] = find(pB{ii});     lin = [lin(:); cB(:)];   end
    if isempty(lin),    continue,   end
    rr = mod(lin-1,nR);     cc = floor((lin-1)/nR);
    occL(mod(rr,NL)+1) = true;  occR(mod(cc,NR)+1) = true;
end
% Bases without the empty degrees; nothing to do when none is empty
[ZLn,ZRn,keepL,keepR] = prune_zero_monomials(ZL,ZR,occL,occR);
if numel(keepL)==NL && numel(keepR)==NR,    return,     end
% Linear positions of the kept entries; A and the columns of B restricted
% to them, in every cell stored at full size
rowK = reshape(keepL(:)+(0:m-1)*NL,[],1);
colK = reshape(keepR(:)+(0:n-1)*NR,[],1);
vecK = reshape(rowK(:)+(colK(:).'-1)*nR,[],1);
for ii = 1:numel(pA)
    if numel(pA{ii})~=nR*nC,    continue,   end
    pA{ii} = pA{ii}(vecK);
    if ~isempty(pB{ii}),    pB{ii} = pB{ii}(:,vecK);   end
end
% The block on the pruned bases
B = sdopvar(struct('A',{pA},'B',{pB}),B.vars,B.Zd,ZLn,ZRn,B.dom,B.dims);
end

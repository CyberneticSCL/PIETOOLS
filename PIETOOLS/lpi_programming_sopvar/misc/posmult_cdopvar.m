function [prog,M,info] = posmult_cdopvar(prog,dims,spaces,dom,w,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,M,INFO] = POSMULT_CDOPVAR(PROG,DIMS,SPACES,DOM,W,OPTIONS) declares a
% positive multiplier on a container space: a self-adjoint 'cdopvar' M that
% acts on every L_2 space as multiplication by the polynomial matrix
%
%    M(theta) = sum_c g_c(theta) Lambda_c(theta)' Q_c Lambda_c(theta),
%               Q_c >= 0,   Lambda_c(theta) = Z_{w_c}(theta) kron I,
%
% one Gram Q_c per Positivstellensatz term c of OPTIONS.psatz, g_c its
% weight (1, the box product, or one face of the box) and w_c the weight
% degree W reduced by OPTIONS.psatz_offset(c). On an R^{m_k} space the block
% is the constant int g_c Lambda_c' Q_c Lambda_c dtheta, and between an R
% and an L_2 space the integral operator with that polynomial kernel: the
% blocks 'copquadvar' builds for its identity and multiplier basis
% operators (alpha = 1), which is what this calls, once per term.
%
% Pointwise, M(theta) >= 0 on the box. In 1-D the pair psatz = [0 1] at
% degrees w and w-1 is the form S0 + (theta-a)(b-theta) S1 of every
% polynomial matrix of degree 2w nonnegative on [a,b] (the interval matrix
% theorem the proof program cites, relevant_positivstellensatz.md Sec. 1;
% scalar case Markov-Lukacs), and the face codes give the odd-degree form
% at degree w. In N-D the terms are the Schmuedgen-type products of face
% weights and no exactness at a fixed degree is claimed.
%
% M is the weight W(theta) of the proof program (proof_roadmap.md eq. (2)),
% declared on its own so that Pop = A'*M*A ('poscopvar_lift') carries the
% certificate as the solved M, readable with 'getsol_lpivar_sop'.
%
% INPUT
% - prog, dims, spaces, dom: as 'poscopvar' (one set of spaces, in = out);
% - w:       weight degree per registry variable: a scalar, a 1 x nv row, or
%            a cell with one such entry per space (per-space degrees give a
%            weight graded by space). Defaults to 1;
% - options: (optional) struct with fields
%   psatz          row of term codes, each 0 (weight 1), 1 (the product
%                  prod_d (theta_d-a_d)(b_d-theta_d)) or 2d+1 / 2d+2 (the
%                  lower / upper face of the box in registry direction d,
%                  normalised, as 'copquadvar'). One Gram per entry;
%                  repeated codes are allowed. Defaults to 0;
%   psatz_offset   degree reduction per term, a scalar or one entry per
%                  code, applied in every direction and floored at 0.
%                  Defaults to 1 for code 1 and 0 otherwise: the degree
%                  count of the interval theorem, deg S1 = deg S0 - 2;
%
% OUTPUT
% - prog:    the program with the Grams declared, in the order of psatz,
%            each with its constraint Q_c >= 0;
% - M:       'cdopvar' on the given spaces, self-adjoint;
% - info:    struct with fields psatz, offset, deg (1 x nterms cell, the
%            per-space 'int' degrees handed to 'poscopvar'), Qcell (1 x
%            nterms cell, as 'poscopvar' returns), gram (1 x nterms Gram
%            dimensions), ndec (decision variables of M).
%
% NOTES
% Cost: one 'poscopvar' call per term (its pair loop over the identity and
% multiplier components, no semiseparable integral) and, for two or more
% terms, their sum by '@cdopvar/plus', an O(q log q) merge of the decision
% variable lists in the q new variables.
% Every block is then pruned to its occupied monomials ('prune_block'):
% copquadvar returns the reachable bases, which pad a multiplier-only
% block's coefficient columns by (2w+2)^nv, measured 64 in 3-D at w = 1.
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
codes = 0;
if isfield(options,'psatz') && ~isempty(options.psatz)
    codes = reshape(double(options.psatz),1,[]);
end
if any(~ismember(codes,[0,1,3:2*nv+2]))
    error("'psatz' entries should be 0, 1, or face codes 2d+1 / 2d+2 with d in 1..nv.")
end
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

M = [];
info = struct('psatz',codes,'offset',off,'deg',{cell(1,nt)},...
              'Qcell',{cell(1,nt)},'gram',zeros(1,nt),'ndec',0);
for c = 1:nt
    deg = cell(1,Ms);
    for k = 1:Ms
        deg{k} = struct('int',max(wsp{k}-off(c),0));
    end
    popt = struct('psatz',codes(c),'include',{incl});
    [prog,Mc,Qc] = poscopvar(prog,dims,spaces,dom,deg,popt);
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
info.ndec = numel(M.Zd);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = expand_row(v,nv)
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
% The block on the monomials that carry a coefficient, per direction, as
% '@sdopvar/mtimes' prunes its result (09/11/2026): a degree is dropped
% only when neither A nor any row of B is nonzero there, so the operator
% is unchanged, and the canonical multiplier form survives (right degree
% 0 stays occupied wherever a multiplier cell has content).
m = B.dims(1);  n = B.dims(2);
ZL = B.ZL;      ZR = B.ZR;
NL = prod([cellfun(@numel,ZL),1]);  NR = prod([cellfun(@numel,ZR),1]);
nR = m*NL;      nC = n*NR;
pA = B.params.A;    pB = B.params.B;
occL = false(NL,1); occR = false(NR,1);
for ii = 1:numel(pA)
    if numel(pA{ii})~=nR*nC,    continue,   end
    lin = find(pA{ii});
    if ~isempty(pB{ii}),    [~,cB] = find(pB{ii});     lin = [lin(:); cB(:)];   end
    if isempty(lin),    continue,   end
    rr = mod(lin-1,nR);     cc = floor((lin-1)/nR);
    occL(mod(rr,NL)+1) = true;  occR(mod(cc,NR)+1) = true;
end
[ZLn,ZRn,keepL,keepR] = prune_zero_monomials(ZL,ZR,occL,occR);
if numel(keepL)==NL && numel(keepR)==NR,    return,     end
rowK = reshape(keepL(:)+(0:m-1)*NL,[],1);
colK = reshape(keepR(:)+(0:n-1)*NR,[],1);
vecK = reshape(rowK(:)+(colK(:).'-1)*nR,[],1);
for ii = 1:numel(pA)
    if numel(pA{ii})~=nR*nC,    continue,   end
    pA{ii} = pA{ii}(vecK);
    if ~isempty(pB{ii}),    pB{ii} = pB{ii}(:,vecK);   end
end
B = sdopvar(struct('A',{pA},'B',{pB}),B.vars,B.Zd,ZLn,ZRn,B.dom,B.dims);
end

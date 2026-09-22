function [opts,report] = eq_opts_sopvar(Dop,opts,ztol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [OPTS,REPORT] = EQ_OPTS_SOPVAR(DOP,OPTS,ZTOL) chooses the 'possopvar'
% options for declaring a positive operator DEOP to enforce DOP + DEOP == 0,
% by inspecting the support of DOP. It is the 'sopvar' counterpart of
% 'get_eq_opts_2D'.
%
% INPUT
% - Dop:   'sopvar' or 'sdopvar' object, the operator being cancelled;
% - opts:  (optional) 'struct' of 'possopvar' options to start from. Any
%          'include' already present is intersected with the surviving
%          blocks, so a caller restriction is never widened;
% - ztol:  (optional) coefficients at or below this magnitude count as zero.
%          Defaults to 1e-12, matching 'get_eq_opts_2D';
%
% OUTPUT
% - opts:   the input options with 'include' set to the basis operators that
%           can contribute and 'sep' set in every direction where DOP is
%           already separable;
% - report: 'struct' with fields 'lines' (cellstr, what was decided and
%           why), 'nkept', 'ndropped' and 'support' (the logical support of
%           DOP over its parameter cells).
%
% WHY BLOCKS CAN BE DROPPED
% 'possopvar' builds DEOP = sum_ij B_i^* Q_ij B_j with Q >= 0. Consider the
% diagonal term B_i^* Q_ii B_i. If every parameter cell it can populate is
% zero in DOP, the equality forces that diagonal block of Q to zero, and a
% positive semidefinite matrix with a zero diagonal block has that block's
% entire row and column zero. Basis operator i therefore contributes
% nothing at all -- not merely nothing useful -- and declaring it only costs
% decision variables and adds identically-zero equality rows.
%
% WHICH CELLS A DIAGONAL REACHES
% MEASURED, by declaring one block at a time and reading the nonzero cells:
% the diagonal of block alpha reaches exactly the product set
%       reach(alpha) = prod_k R_k,   R_k = {1} if alpha_k==1 else {2,3},
% because adjoining maps a lower integral to an upper one, so a direction
% carrying an integral contributes to BOTH Volterra cells while a multiplier
% direction contributes only to the delta cell. A block is dropped only when
% DOP is zero on all of reach(alpha).
%
% This is deliberately weaker than 'get_eq_opts_2D', which drops component
% (i,j) when cell (i,j) alone is zero. That is not implied by the argument
% above: the diagonal of block (2,2) also reaches cells (3,2), (2,3) and
% (3,3), so its Gram block is forced to zero only when all four vanish.
% Dropping a block whose diagonal can still reach a nonzero cell restricts
% the cone, which cannot give a wrong answer but can lose a certificate.
%
% SEPARABILITY
% Where DOP already has equal lower and upper kernels in a direction -- the
% same test 'get_eq_opts_2D' applies -- that direction is marked separable,
% so 'possopvar' uses one full-domain integral in place of two Volterra
% ones. That halves the blocks in the direction at no cost in reach.
%
% See also POSSOPVAR, GET_EQ_OPTS_2D, SETTINGS2POSSOPVAR, LPI_EQ_SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/21/2026: the negativity operator dominates the cost
%                  of an LPI -- at two spatial variables its Gram is the
%                  overwhelming majority of the program -- and most of its
%                  basis is dead on arrival, because the operator being
%                  cancelled has no content in the cells those blocks can
%                  reach. Declaring them costs decision variables and
%                  contributes only identically-zero equality rows. The
%                  'opvar2d' path has 'get_eq_opts_2D' for this; 'sopvar'
%                  had nothing.

if nargin<2 || isempty(opts)
    opts = struct();
end
if nargin<3 || isempty(ztol)
    ztol = 1e-12;
end
if ~isa(Dop,'sopvar') && ~isa(Dop,'sdopvar')
    error("The operator should be an 'sopvar' or 'sdopvar' object.")
end

report = struct('lines',{{}},'nkept',0,'ndropped',0,'support',[]);

% % % Support of Dop over its parameter cells. A cell of an 'sdopvar' counts
% as zero only when neither the constant part nor any row of the decision
% part has an entry above the tolerance, i.e. when the coefficient vanishes
% for EVERY value of the decision variables.
ncell = numcells(Dop);
n3 = round(log(ncell)/log(3));
if 3^n3 ~= ncell
    error("The parameter cell count %d is not a power of 3.",ncell)
end
sup = false(1,ncell);
for k = 1:ncell
    sup(k) = cellmax(Dop,k) > ztol;
end
report.support = sup;
if ~any(sup)
    report.lines{end+1} = ['Dop is identically zero; every basis operator ' ...
        'is dropped except the multiplier, which is kept so the result ' ...
        'remains a declarable operator.'];
end

% % % Directions in which Dop is already separable: the gamma_k=2 and
% gamma_k=3 cells agree for every setting of the other indices. Same test as
% 'get_eq_opts_2D', which reads it off Qop.R22.
sepv = false(1,n3);
if isfield(opts,'sep') && ~isempty(opts.sep)
    sepv = logical(reshape(opts.sep,1,[]));
    if isscalar(sepv), sepv = repmat(sepv,1,n3); end
end
for kdir = 1:n3
    if sepv(kdir), continue, end
    ok = true;
    for k = 1:ncell
        g = gamma_of(k,n3);
        if g(kdir)~=2, continue, end
        g3 = g;  g3(kdir) = 3;
        if celldiff(Dop,k,lin_of(g3,n3)) > ztol
            ok = false;  break
        end
    end
    if ok && any(sup)
        sepv(kdir) = true;
        report.lines{end+1} = sprintf(['direction %d: Dop has equal lower ' ...
            'and upper kernels, so it is declared separable -- one ' ...
            'full-domain integral in place of two Volterra ones.'],kdir);
    end
end

% % % Candidate blocks, and the cells each diagonal can reach.
alpha_all = enum_alpha(n3,sepv);
keep = false(size(alpha_all,1),1);
for r = 1:size(alpha_all,1)
    cells = reach(alpha_all(r,:),n3);
    keep(r) = any(sup(cells));
end

% A degenerate all-zero target would leave nothing declarable; keep the
% all-multiplier block so the caller still gets a valid operator.
if ~any(keep)
    keep(1) = true;
end

incl = alpha_all(keep,:);
if isfield(opts,'include') && ~isempty(opts.include) ...
        && size(opts.include,2)==n3
    % Never widen a restriction the caller already asked for.
    incl = incl(ismember(incl,opts.include,'rows'),:);
    if isempty(incl), incl = opts.include(1,:); end
end

opts.include = incl;
opts.sep = sepv;
report.nkept = size(incl,1);
report.ndropped = size(alpha_all,1) - size(incl,1);
report.lines{end+1} = sprintf(['kept %d of %d basis operators; %d dropped ' ...
    'because Dop is zero on every cell their diagonal can reach.'], ...
    report.nkept,size(alpha_all,1),report.ndropped);

end


%%
function C = reach(a,n3)
% Linear cell indices the diagonal of block 'a' can populate:
%   reach(a) = prod_k ({1} if a_k==1 else {2,3}).
% MEASURED by declaring one block at a time; adjoining sends a lower
% integral to an upper one, so an integral direction reaches both Volterra
% cells and a multiplier direction only the delta cell. alpha_k = 4, the
% full-domain integral, is an integral in this respect.
sets = cell(1,n3);
for k = 1:n3
    if a(k)==1
        sets{k} = 1;
    else
        sets{k} = [2 3];
    end
end
G = sets{1}(:);
for k = 2:n3
    G = [kron(ones(numel(sets{k}),1),G), ...
         kron(sets{k}(:),ones(size(G,1),1))];
end
C = zeros(size(G,1),1);
for r = 1:size(G,1)
    C(r) = lin_of(G(r,:),n3);
end
end


%%
function A = enum_alpha(n3,sepv)
% Admissible multi-indices: {1,2,3} in an ordinary direction and {1,4} in a
% separable one, direction 1 varying fastest, matching 'possopvar'.
vals = cell(1,n3);
for k = 1:n3
    if sepv(k), vals{k} = [1,4]; else, vals{k} = [1,2,3]; end
end
A = vals{1}(:);
for k = 2:n3
    A = [kron(ones(numel(vals{k}),1),A), ...
         kron(vals{k}(:),ones(size(A,1),1))];
end
end


%%
function n = numcells(P)
if isa(P,'sdopvar'), n = numel(P.params.A); else, n = numel(P.params); end
end


%%
function m = cellmax(P,k)
% Largest magnitude in cell k, over the constant part and every row of the
% decision part.
m = 0;
if isa(P,'sdopvar')
    if k<=numel(P.params.A)
        m = max(m,absmax(P.params.A{k}));
        m = max(m,absmax(P.params.B{k}));
    end
else
    if k<=numel(P.params)
        m = absmax(P.params{k});
    end
end
end


%%
function d = celldiff(P,k1,k2)
% Largest magnitude of the difference between two cells. Cells of an
% 'sopvar' share one global ZL/ZR, so they are directly comparable.
if isa(P,'sdopvar')
    d = max(absmax(getA(P,k1)-getA(P,k2)), absmax(getB(P,k1)-getB(P,k2)));
else
    d = absmax(getP(P,k1)-getP(P,k2));
end
end

function X = getA(P,k)
if k<=numel(P.params.A), X = P.params.A{k}; else, X = 0; end
end
function X = getB(P,k)
if k<=numel(P.params.B), X = P.params.B{k}; else, X = 0; end
end
function X = getP(P,k)
if k<=numel(P.params), X = P.params{k}; else, X = 0; end
end

function m = absmax(X)
if isempty(X), m = 0; else, m = full(max(abs([nonzeros(X);0]))); end
end


%%
function g = gamma_of(k,n3)
g = ones(1,n3);
for t = 1:n3
    g(t) = mod(floor((k-1)/3^(t-1)),3)+1;
end
end

function k = lin_of(g,n3)
k = 1;
for t = 1:n3
    k = k + (g(t)-1)*3^(t-1);
end
end

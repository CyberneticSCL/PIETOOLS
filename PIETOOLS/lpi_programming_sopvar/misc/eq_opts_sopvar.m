function [opts,report] = eq_opts_sopvar(Dop,opts,ztol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [OPTS,REPORT] = EQ_OPTS_SOPVAR(DOP,OPTS,ZTOL) sets the options sep and
% include of POSSOPVAR for a positive operator Deop declared to satisfy
% DOP + Deop = 0, from the support of DOP. With n3 the number of variables
% common to the input and output of DOP, gamma over its 3^n3 parameter
% cells (gamma_k = 1 multiplier, 2 lower integral, 3 upper integral), and
% a cell nonzero when an entry of its coefficients (for an 'sdopvar', of A
% or of a row of B) exceeds ZTOL in magnitude:
%   - direction k is separable when DOP has a nonzero cell and, for every
%     gamma with gamma_k = 2, the cells gamma and gamma + e_k differ by at
%     most ZTOL in every entry;
%   - a basis operator alpha, alpha_k in {1,2,3}, or {1,4} in a separable
%     direction, is kept when DOP is nonzero on a cell of
%
%         reach(alpha) = R_1 x ... x R_n3,
%         R_k = {1} for alpha_k = 1,   R_k = {2,3} otherwise,
%
%     and alpha = (1,...,1) is kept when none is.
% A direction with OPTS.sep true is separable without the test, and an
% OPTS.include of rows alpha restricts the kept rows to its own (to its
% first row when none remains).
%
% INPUT
% - Dop:     'sopvar' or 'sdopvar' object with 3^n3 parameter cells;
% - opts:    (optional) 'struct' of POSSOPVAR options (default struct());
%            read are sep ('logical', scalar or 1 x n3) and include (a
%            'double' array of rows alpha with n3 columns; an include of
%            another form is replaced);
% - ztol:    (optional) 'double', the largest magnitude counted as zero
%            (default 1e-12).
%
% OUTPUT
% - opts:    'struct', OPTS with the fields include ('double', the kept
%            rows alpha, n3 columns, in the order of the enumeration of
%            every alpha, direction 1 fastest) and sep (1 x n3 'logical',
%            the separable directions);
% - report:  'struct' with fields lines ('cellstr', one sentence per
%            decision), nkept and ndropped ('double', the numbers of kept
%            and dropped basis operators) and support (1 x 3^n3 'logical',
%            true on the nonzero cells of DOP, in linear cell order,
%            direction 1 fastest).
%
% See also POSSOPVAR, GET_EQ_OPTS_2D, SETTINGS2POSSOPVAR, LPI_EQ_SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eq_opts_sopvar
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% MMP, 09/30/2026: The parameter-cell <-> gamma map (Sec. 1 and 4 of
%                  sopvar_implementation_notes.pdf; direction 1 fastest)
%                  and the multi-index enumerations
%                  call the shared 'gamma_of_cell', 'cell_of_gamma' and
%                  'multiindex_grid' (sopvar/misc/conventions); the locals
%                  gamma_of and lin_of, and the kron loops of reach and
%                  enum_alpha, are gone. The cell order is load-bearing: the
%                  'include' rows emitted here are read by 'copquadvar' and
%                  'sopquadvar', so all three must enumerate alike; the
%                  helper states that order once.
%                  Same output on every input that ran before; at n3 = 0
%                  (no shared variable) reach and enum_alpha now return
%                  cell 1 and the empty multi-index where they errored.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), the name line of the licence block
%                  added, step comments in the body and headers of the
%                  local functions; code unchanged. The deleted header text
%                  was the argument that a dropped block contributes
%                  nothing, how the reach sets were found (by measurement)
%                  and a comparison with 'get_eq_opts_2D'.

% Defaults: OPTS struct(), ZTOL 1e-12; DOP an 'sopvar' or 'sdopvar'
if nargin<2 || isempty(opts)
    opts = struct();
end
if nargin<3 || isempty(ztol)
    ztol = 1e-12;
end
if ~isa(Dop,'sopvar') && ~isa(Dop,'sdopvar')
    error("The operator should be an 'sopvar' or 'sdopvar' object.")
end

% Empty report
report = struct('lines',{{}},'nkept',0,'ndropped',0,'support',[]);

% % % Support of Dop over its parameter cells. A cell of an 'sdopvar' counts
% as zero only when neither the constant part nor any row of the decision
% part has an entry above the tolerance, i.e. when the coefficient vanishes
% for EVERY value of the decision variables.
% n3 from the cell count 3^n3; then the support, cell by cell
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
% A zero DOP, reported
if ~any(sup)
    report.lines{end+1} = ['Dop is identically zero; every basis operator ' ...
        'is dropped except the multiplier, which is kept so the result ' ...
        'remains a declarable operator.'];
end

% % % Directions in which Dop is already separable: the gamma_k=2 and
% gamma_k=3 cells agree for every setting of the other indices. Same test as
% 'get_eq_opts_2D', which reads it off Qop.R22.
% The directions OPTS.sep marks are separable without the test
sepv = false(1,n3);
if isfield(opts,'sep') && ~isempty(opts.sep)
    sepv = logical(reshape(opts.sep,1,[]));
    if isscalar(sepv), sepv = repmat(sepv,1,n3); end
end
% Every other direction: each gamma_k = 2 cell against its gamma_k = 3
% partner; separable when all agree and DOP is nonzero
for kdir = 1:n3
    if sepv(kdir), continue, end
    ok = true;
    for k = 1:ncell
%       g = gamma_of(k,n3);                                                 % MMP, 09/30/2026 (was)
        g = gamma_of_cell(k,n3);                                            % MMP, 09/30/2026
        if g(kdir)~=2, continue, end
        g3 = g;  g3(kdir) = 3;
%       if celldiff(Dop,k,lin_of(g3,n3)) > ztol                             % MMP, 09/30/2026 (was)
        if celldiff(Dop,k,cell_of_gamma(g3)) > ztol                         % MMP, 09/30/2026
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

% The kept rows, within the rows of a caller OPTS.include
incl = alpha_all(keep,:);
if isfield(opts,'include') && ~isempty(opts.include) ...
        && size(opts.include,2)==n3
    % Never widen a restriction the caller already asked for.
    incl = incl(ismember(incl,opts.include,'rows'),:);
    if isempty(incl), incl = opts.include(1,:); end
end

% The options and the counts of the report
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
% INPUT  a:  1 x n3 'double', the multi-index alpha;
%        n3: 'double', the number of variables.
% OUTPUT C:  'double' column of linear cell indices.
% Per direction {1} or {2,3}
sets = cell(1,n3);
for k = 1:n3
    if a(k)==1
        sets{k} = 1;
    else
        sets{k} = [2 3];
    end
end
% G = sets{1}(:);                                                           % MMP, 09/30/2026 (was)
% for k = 2:n3                                                              % MMP, 09/30/2026 (was)
%     G = [kron(ones(numel(sets{k}),1),G), ...
%          kron(sets{k}(:),ones(size(G,1),1))];                             % MMP, 09/30/2026 (was)
% end                                                                       % MMP, 09/30/2026 (was)
% C = zeros(size(G,1),1);                                                   % MMP, 09/30/2026 (was)
% for r = 1:size(G,1)                                                       % MMP, 09/30/2026 (was)
%     C(r) = lin_of(G(r,:),n3);                                             % MMP, 09/30/2026 (was)
% end                                                                       % MMP, 09/30/2026 (was)
% The product set, direction 1 fastest, as linear cell indices (column).    % MMP, 09/30/2026
C = cell_of_gamma(multiindex_grid(sets));                                   % MMP, 09/30/2026
end


%%
function A = enum_alpha(n3,sepv)
% Admissible multi-indices: {1,2,3} in an ordinary direction and {1,4} in a
% separable one, direction 1 varying fastest, matching 'possopvar'.
% INPUT  n3:   'double', the number of variables;
%        sepv: 1 x n3 'logical', the separable directions.
% OUTPUT A:    'double', one multi-index per row, n3 columns.
% Admissible entries per direction, then their grid
vals = cell(1,n3);
for k = 1:n3
    if sepv(k), vals{k} = [1,4]; else, vals{k} = [1,2,3]; end
end
% A = vals{1}(:);                                                           % MMP, 09/30/2026 (was)
% for k = 2:n3                                                              % MMP, 09/30/2026 (was)
%     A = [kron(ones(numel(vals{k}),1),A), ...
%          kron(vals{k}(:),ones(size(A,1),1))];                             % MMP, 09/30/2026 (was)
% end                                                                       % MMP, 09/30/2026 (was)
A = multiindex_grid(vals);                                                  % MMP, 09/30/2026
end


%%
function n = numcells(P)
% N = NUMCELLS(P) returns the number of parameter cells of P.
% INPUT  P: 'sopvar' or 'sdopvar' object.
% OUTPUT n: 'double'.
if isa(P,'sdopvar'), n = numel(P.params.A); else, n = numel(P.params); end
end


%%
function m = cellmax(P,k)
% Largest magnitude in cell k, over the constant part and every row of the
% decision part.
% INPUT  P: 'sopvar' or 'sdopvar' object; k: 'double', a cell index.
% OUTPUT m: 'double', 0 for a cell beyond the last.
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
% INPUT  P:      'sopvar' or 'sdopvar' object;
%        k1, k2: 'double', two cell indices.
% OUTPUT d:      'double'; for an 'sdopvar', the larger of the A and B
%                differences.
if isa(P,'sdopvar')
    d = max(absmax(getA(P,k1)-getA(P,k2)), absmax(getB(P,k1)-getB(P,k2)));
else
    d = absmax(getP(P,k1)-getP(P,k2));
end
end

function X = getA(P,k)
% X = GETA(P,K) returns the constant part A of cell K of an 'sdopvar', 0
% beyond the last cell.
% INPUT  P: 'sdopvar' object; k: 'double', a cell index.
% OUTPUT X: 'double' matrix, or 0.
if k<=numel(P.params.A), X = P.params.A{k}; else, X = 0; end
end
function X = getB(P,k)
% X = GETB(P,K) returns the decision part B of cell K of an 'sdopvar', 0
% beyond the last cell.
% INPUT  P: 'sdopvar' object; k: 'double', a cell index.
% OUTPUT X: 'double' matrix, or 0.
if k<=numel(P.params.B), X = P.params.B{k}; else, X = 0; end
end
function X = getP(P,k)
% X = GETP(P,K) returns cell K of an 'sopvar', 0 beyond the last cell.
% INPUT  P: 'sopvar' object; k: 'double', a cell index.
% OUTPUT X: 'double' matrix, or 0.
if k<=numel(P.params), X = P.params{k}; else, X = 0; end
end

function m = absmax(X)
% M = ABSMAX(X) returns the largest magnitude of the entries of X, 0 for an
% empty or zero X.
% INPUT  X: 'double' matrix, full or sparse.
% OUTPUT m: 'double', full.
if isempty(X), m = 0; else, m = full(max(abs([nonzeros(X);0]))); end
end


% MMP, 09/30/2026: local functions gamma_of and lin_of (Initial coding
% 09/21/2026, no line stamps) deleted: 'gamma_of_cell' and 'cell_of_gamma'
% (sopvar/misc/conventions) compute the same map.

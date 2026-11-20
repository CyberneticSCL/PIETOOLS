function [deg,report] = degbalance_core(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,REPORT] = DEGBALANCE_CORE(P,OPTS) picks the degrees of the basis of
% a positive operator so that, declared with them, it has monomials of
% similar degree to those of the operator P.
%
% INPUT
% - P:    'sopvar' or 'sdopvar' object, the operator whose degrees are to be
%         matched -- typically the one a positive operator must cancel;
% - opts: (optional) 'struct' with fields
%         - include: nblk x n3 array, the basis multi-indices to size, one
%                    row per block, an entry 1, 2 or 3 per direction (1 or 4
%                    in a separable direction). Defaults to every block, or,
%                    when 'sep' is given, every block over {1,4} per
%                    separable direction;
%         - sep:     1 x n3 logical, the separable directions;
%         - weight:  'uniform' (default) halves every cap. 'slotwise'
%                    instead halves only the integration slot and leaves the
%                    output slot at full degree -- see SLOT WEIGHTING;
%
% OUTPUT
% - deg:    1 x nblk cell, one degree struct per row of opts.include, with
%           fields 'int' (1 x n3, the integration-variable caps), 'mult'
%           (1 x n3, the output-variable caps), 'joint' (the cap of the full
%           set of variables) and 'subset' (1 x 2^(2*n3), the cap of every
%           subset of variables);
% - report: 'struct' with 'lines' (what was decided) and 'maxdeg' (the
%           per-block per-subset maximal degrees read off P).
%
% THE RULE
% For every block, halve the maximal degree of the target on every
% nonempty subset of variables,
%       d{i,j} = ceil(maxdeg{i,j}/2),
% the degrees read from the coefficient support of P together with ZL and
% ZR. The halving is the Newton-polytope factor: a positive operator is a
% sum of squares in the Gram sense, so its basis carries roughly half the
% degree of what it represents.
%
% SLOT WEIGHTING
% The two slots do not contribute equally. MEASURED on a one-variable
% Volterra block, the maximal degree of the DECLARED operator is
%       2*deg.int + deg.mult + 1,
% because the integrand carries r^(a_p+c_q) from both factors, so the
% variable-limit integral int_a^s gives s^(a_p+c_q+1), on top of the
% s^(b_p) of the left factor. The integration slot therefore reaches twice
% as far per unit of degree as the output slot. 'uniform' ignores this and
% is conservative on the output slot; 'slotwise' halves only the
% integration slot, as the measurement supports; it has not been validated
% on a solve.
%
% LIMITS
% It sizes the basis; it does not prune it, nor decide which blocks can
% contribute. Per-monomial pruning does not follow from the same argument,
% so the degrees here are a matching rule, not a minimal one.
%
% See also DEGBALANCE, EQ_OPTS_SOPVAR, POSSOPVAR, SETTINGS2POSSOPVAR.
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
% Initial coding MMP, 09/21/2026

% Defaults and checks of the inputs; an operator with no shared variable
% (n3 = 0) is not supported
if nargin<2 || isempty(opts)
    opts = struct();
end
if ~isa(P,'sopvar') && ~isa(P,'sdopvar')
    error("The operator should be an 'sopvar' or 'sdopvar' object.")
end
weight = getdef(opts,'weight','uniform');
if ~ismember(weight,{'uniform','slotwise'})
    error("'weight' should be 'uniform' or 'slotwise'.")
end

% The number of parameter cells of the block P, and the basis sizes
if isa(P,'sdopvar'), ncell = numel(P.params.A); else, ncell = numel(P.params); end
n3 = round(log(ncell)/log(3));
if 3^n3 ~= ncell
    error("The parameter cell count %d is not a power of 3.",ncell)
end
nL = cellfun(@numel,P.ZL);      nL = nL(:).';
nR = cellfun(@numel,P.ZR);      nR = nR(:).';
NL = prod([nL,1]);              NR = prod([nR,1]);

% Separable directions, one flag per shared variable
sepv = getdef(opts,'sep',false(1,n3));
sepv = logical(reshape(sepv,1,[]));
if isscalar(sepv), sepv = repmat(sepv,1,n3); end

incl = getdef(opts,'include',[]);
if isempty(incl)
    % every basis multi-index, {1,2,3} per direction or {1,4} in a
    % separable one, direction 1 fastest
    vals = cell(1,n3);
    for k = 1:n3
        if sepv(k), vals{k} = [1,4]; else, vals{k} = [1,2,3]; end
    end
    incl = multiindex_grid(vals);          % direction 1 fastest
end

% One degree struct per block
report = struct('lines',{{}},'maxdeg',{cell(1,size(incl,1))});
deg = cell(1,size(incl,1));

% Degree tables of the two bases: row t of DL gives the degree of composite
% ZL monomial t in each output variable, first direction outermost.
DL = multiindex_grid(P.ZL,'first_slowest');
DR = multiindex_grid(P.ZR,'first_slowest');

for r = 1:size(incl,1)
    a = incl(r,:);
    % The cells the diagonal of block a can populate: every gamma in the
    % product of the sets {1} or {2,3} per direction, direction 1 fastest,
    % as cells
    sets = cell(1,n3);
    for k = 1:n3
        if a(k)==1, sets{k} = 1; else, sets{k} = [2 3]; end
    end
    cells = cell_of_gamma(multiindex_grid(sets));

    % % % Exponent pairs the target actually carries on those cells. The
    % rows of a coefficient are (matrix row outer, ZL monomial inner) and
    % the columns likewise over ZR.
    EO = zeros(0,n3);   EI = zeros(0,n3);
    for c = cells(:).'
        [ri,ci] = nzpos(P,c);
        if isempty(ri), continue, end
        EO = [EO; DL(mod(ri-1,NL)+1,:)];                             %#ok<AGROW>
        EI = [EI; DR(mod(ci-1,NR)+1,:)];                             %#ok<AGROW>
    end

    % % % Per-subset maximal degrees of the target, over the variable order
    % [out_1..out_n3, in_1..in_n3]; subset i is the set of bits of i-1.
    nv = 2*n3;
    T = zeros(1,2^nv);
    if ~isempty(EO)
        E = [EO, EI];
        for i = 2:2^nv
            b = bitget(i-1,1:nv)>0;
            T(i) = max(sum(E(:,b),2));
        end
    end
    report.maxdeg{r} = T;

    % % % Halve. 'uniform' halves every entry; 'slotwise' halves only the
    % subsets made purely of integration slots and leaves the rest whole,
    % reflecting that the integration slot reaches twice as far.
    % CEIL, not floor: under-sizing the basis shows up as spurious
    % infeasibility -- the positive operator cannot represent the target, so
    % the equality forces genuine coefficients to zero -- while over-sizing
    % only costs decision variables.
    S = zeros(1,2^nv);
    for i = 2:2^nv
        b = bitget(i-1,1:nv)>0;
        if strcmp(weight,'uniform') || all(b(1:n3) | ~b)
            S(i) = ceil(T(i)/2);
        else
            S(i) = T(i);
        end
    end

    % Singletons and the full set are the int, mult and joint caps. SLOT
    % CORRESPONDENCE: 'int', the INTEGRATION-variable cap of the basis, takes
    % the degree of the INPUT variable of the target, and 'mult', the
    % output-variable cap, the degree of its OUTPUT variable.
    d = struct();
    d.int  = zeros(1,n3);
    d.mult = zeros(1,n3);
    for k = 1:n3
        d.int(k)  = S(1+2^(n3+k-1));    % from the target's input slot
        d.mult(k) = S(1+2^(k-1));       % from the target's output slot
    end
    if a(1)~=0
        d.mult(a==1) = 0;   % a multiplier direction carries no output slot
    end
    d.joint  = S(2^nv);
    d.subset = S;
    deg{r} = d;
end

report.lines{end+1} = sprintf(['sized %d blocks from the target degrees, ' ...
    'weight ''%s''.'],size(incl,1),weight);

end


%%
function [ri,ci] = nzpos(P,k)
% [RI,CI] = NZPOS(P,K) returns the row and column positions of the stored
% nonzeros of cell K of P (constant or decision part), empty past the last
% cell.
% INPUT   P: 'sopvar' or 'sdopvar'; k: linear cell index.
% OUTPUT  ri, ci: column vectors of positions in the unvectorized C_k.
if isa(P,'sdopvar')
    if k>numel(P.params.A), ri=[]; ci=[]; return, end
    m = size(P.params.A{k},1);
    [r1,~] = find(P.params.A{k});
    [~,c2] = find(P.params.B{k});
    lin = unique([r1(:);c2(:)]);
    nrow = P.dims(1)*prod([cellfun(@numel,P.ZL),1]);
    ri = mod(lin-1,nrow)+1;
    ci = floor((lin-1)/nrow)+1;
    if m==0, ri=[]; ci=[]; end
else
    if k>numel(P.params), ri=[]; ci=[]; return, end
    [ri,ci] = find(P.params{k});
end
end


%%
function v = getdef(s,f,dflt)
% V = GETDEF(S,F,DFLT) reads field F of the struct S, or DFLT when it is
% absent or empty.
if isa(s,'struct') && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end

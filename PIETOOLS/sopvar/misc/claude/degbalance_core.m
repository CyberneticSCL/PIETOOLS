function [deg,report] = degbalance_core(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,REPORT] = DEGBALANCE_CORE(P,OPTS) picks 'possopvar' degrees so that a
% positive operator declared with them has monomials of similar degree to
% those of P. It is the n-variate counterpart of '@opvar/degbalance' and
% '@dopvar2d/degbalance', and is shared by '@sopvar/degbalance' and
% '@sdopvar/degbalance'.
%
% INPUT
% - P:    'sopvar' or 'sdopvar' object, the operator whose degrees are to be
%         matched -- typically the one a positive operator must cancel;
% - opts: (optional) 'struct' with fields
%         - include: the basis multi-indices to size, as 'possopvar' takes
%                    them. Defaults to every block, or, when 'sep' is given,
%                    every block over {1,4} per separable direction. Pass
%                    the 'include' from 'eq_opts_sopvar' to size only the
%                    blocks that can contribute;
%         - sep:     1 x n3 logical, the separable directions;
%         - weight:  'uniform' (default) halves every cap, which is what
%                    '@dopvar2d/degbalance' does as shipped. 'slotwise'
%                    instead halves only the integration slot and leaves the
%                    output slot at full degree -- see the note below;
%
% OUTPUT
% - deg:    1 x nblk cell, one degree struct per row of opts.include, with
%           fields 'int', 'mult', 'joint' and 'subset', ready to hand to
%           'possopvar' as its degree argument;
% - report: 'struct' with 'lines' (what was decided) and 'maxdeg' (the
%           per-block per-subset maximal degrees read off P).
%
% THE RULE
% '@dopvar2d/degbalance' as shipped is, for every block,
%       d2{i,j} = floor(Pmaxdegs.R22{i,j}/2),
% i.e. halve the target's per-NONEMPTY-SUBSET maximal degree array
% elementwise. That form is already n-variate, since the array is indexed by
% variable subset, so nothing needs generalizing beyond reading the array off
% an 'sopvar' -- which has no 'getdeg', so the degrees are taken here from
% the coefficient support together with ZL and ZR.
%
% The halving is the Newton-polytope factor: a positive operator is a sum of
% squares in the Gram sense, so its basis carries roughly half the degree of
% what it represents.
%
% SLOT WEIGHTING
% The two slots do not contribute equally. MEASURED on a one-variable
% Volterra block, the maximal degree of the DECLARED operator is
%       2*deg.int + deg.mult + 1,
% because the integrand carries r^(a_p+c_q) from both factors, so the
% variable-limit integral int_a^s gives s^(a_p+c_q+1), on top of the left
% factor's s^(b_p). The integration slot therefore reaches twice as far per
% unit of degree as the output slot. 'uniform' ignores this and is
% conservative on the output slot; 'slotwise' halves only the integration
% slot. '@dopvar2d/degbalance' carries the same observation in a
% commented-out alternative ("degrees in the primary variables get doubled,
% whereas those in the dummy variables do not") with factors 0.5 and 1,
% which agree with the measurement. 'uniform' is the default because it is
% what the toolbox ships and has been exercised; 'slotwise' is offered
% because the measurement supports it, not because it has been validated on
% a solve.
%
% WHAT THIS DOES NOT DO
% It sizes the basis; it does not prune it. Deciding WHICH blocks can
% contribute is 'eq_opts_sopvar'. Per-monomial pruning does not follow from
% the same argument -- see the note in 'eq_opts_sopvar' -- so the degrees
% here are a matching rule, not a minimal one.
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
% Initial coding MMP, 09/21/2026: 'degbalance' existed for 'opvar',
%                  'dopvar' and 'dopvar2d' but not for 'sopvar', so a
%                  positive operator for an 'sopvar' equality had to be
%                  given degrees by hand. The shipped 2D rule is already
%                  n-variate in form -- halve the per-subset maximal degree
%                  array -- so what was missing was reading that array off
%                  an 'sopvar', which has no 'getdeg'.

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

ncell = numcells(P);
n3 = round(log(ncell)/log(3));
if 3^n3 ~= ncell
    error("The parameter cell count %d is not a power of 3.",ncell)
end
nL = cellfun(@numel,P.ZL);      nL = nL(:).';
nR = cellfun(@numel,P.ZR);      nR = nR(:).';
NL = prod([nL,1]);              NR = prod([nR,1]);

sepv = getdef(opts,'sep',false(1,n3));
sepv = logical(reshape(sepv,1,[]));
if isscalar(sepv), sepv = repmat(sepv,1,n3); end

incl = getdef(opts,'include',[]);
if isempty(incl)
    incl = enum_alpha(n3,sepv);
end

report = struct('lines',{{}},'maxdeg',{cell(1,size(incl,1))});
deg = cell(1,size(incl,1));

% Degree tables of the two bases: row t of DL gives the degree of composite
% ZL monomial t in each output variable, first direction outermost.
DL = degree_table(P.ZL);
DR = degree_table(P.ZR);

for r = 1:size(incl,1)
    a = incl(r,:);
    cells = reach(a,n3);

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
    % [out_1..out_n3, in_1..in_n3], which is the order the basis uses once
    % the target's output variable is matched to the basis integration slot
    % and its input variable to the basis output slot. That is the positional
    % correspondence '@dopvar2d/degbalance' uses, and the covering test in
    % the verification script is what confirms it.
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
    % CEIL, not floor. The two shipped implementations disagree:
    % '@opvar/degbalance' rounds up and '@dopvar2d/degbalance' rounds down.
    % Rounding up is taken here because under-sizing the basis shows up as
    % spurious infeasibility -- the positive operator cannot represent the
    % target, so the equality forces genuine coefficients to zero -- while
    % over-sizing only costs decision variables. It is also what makes the
    % cross-check against '@opvar/degbalance' agree exactly.
    S = zeros(1,2^nv);
    for i = 2:2^nv
        b = bitget(i-1,1:nv)>0;
        if strcmp(weight,'uniform') || all(b(1:n3) | ~b)
            S(i) = ceil(T(i)/2);
        else
            S(i) = T(i);
        end
    end

    % Singletons and the full set are the int, mult and joint caps.
    %
    % SLOT CORRESPONDENCE, determined by cross-check against
    % '@opvar/degbalance' rather than from documentation, which is wrong on
    % this point in 'possopvar' itself. Their d{2} is ordered
    % [deg s, deg theta, joint]; matching their answer on targets with
    % DIFFERENT s and theta degrees shows that possopvar's 'int' -- the
    % INTEGRATION-variable cap -- takes the target's INPUT-variable degree,
    % and 'mult' -- the output-variable cap -- takes the target's OUTPUT
    % degree. Reading it the other way round reproduced their numbers
    % transposed on 3 of 4 test targets.
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
function D = degree_table(Z)
% One row per composite monomial of kron(Z{1},...,Z{N}), one column per
% direction, FIRST direction outermost -- the convention 'monomial_gather'
% documents and the class stores.
D = zeros(1,0);
for k = 1:numel(Z)
    D = [kron(D,ones(numel(Z{k}),1)), ...
         kron(ones(max(size(D,1),1),1),Z{k}(:))];
end
end


%%
function C = reach(a,n3)
% Cells the diagonal of block 'a' can populate; see 'eq_opts_sopvar', where
% this was measured.
sets = cell(1,n3);
for k = 1:n3
    if a(k)==1, sets{k} = 1; else, sets{k} = [2 3]; end
end
G = sets{1}(:);
for k = 2:n3
    G = [kron(ones(numel(sets{k}),1),G), ...
         kron(sets{k}(:),ones(size(G,1),1))];
end
C = zeros(size(G,1),1);
for r = 1:size(G,1)
    k = 1;
    for t = 1:n3, k = k + (G(r,t)-1)*3^(t-1); end
    C(r) = k;
end
end


%%
function A = enum_alpha(n3,sepv)
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
function [ri,ci] = nzpos(P,k)
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
function n = numcells(P)
if isa(P,'sdopvar'), n = numel(P.params.A); else, n = numel(P.params); end
end


%%
function v = getdef(s,f,dflt)
if isa(s,'struct') && isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = dflt;
end
end

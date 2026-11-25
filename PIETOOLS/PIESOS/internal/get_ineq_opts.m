function [Z,opts,kst] = get_ineq_opts(f,opts)
% [Z,OPTS,KST] = GET_INEQ_OPTS(F,OPTS) returns the distributed monomial
% basis Z and the options of `piesos_sosvar` for the SOS functional
%   g = <Zop*Z,Qmat*Zop*Z>
% by which `piesos_ineq` imposes F>=0: the type of output of Zop for
% opts.otype='auto', and the monomial degrees and term types of Zop, set
% from the kernels of F unless specified in OPTS.
%
% INPUTS
% - f:      'polyopvar' object, a scalar functional (f.matdim=[1,1], no
%           output variables) with 'intop' coefficients, or constants for
%           degree-0 monomials, cleaned by `clean_poly`, as in
%           `piesos_ineq`;
% - opts:   (optional) struct, the options of `piesos_ineq`. Used here:
%           otype, degs (or deg), exclude, psatz;
%
% OUTPUTS
% - Z:      nZ x 1 'polyopvar' object, the distributed monomial basis of g:
%           the monomials of half the degrees of those of f, within half
%           the Newton polytope of f, without monomials mixing states that
%           depend on different spatial variables (discarded with a
%           warning, id PIESOS:ineq:mixedDependence). With tensor-product
%           rows (otype 2 or 3), each monomial of different states also in
%           each other order of its factors (dmonomials with 'full'). [] if
%           no monomial remains;
% - opts:   struct, OPTS with the fields
%    otype:   1, 2 or 3 if OPTS.otype was 'auto', see NOTES;
%    degs:    if neither degs nor deg is specified: struct with fields 'in'
%             (nf x N array, the degree d of the dummy variables of Zop in
%             each of its nf factors) and 'out' (1 x N array, the degree p
%             of the primary variables in Qmat), see the NOTES of
%             `piesos_ineq`; 0 if Z has no spatial variables;
%    exclude: if not specified and Z has spatial variables: 1 x q cell,
%             the option exclude of `monopmat` for each block of Zop (the
%             term types not needed per variable);
% - kst:    struct, the statistics of the exponents of the monomials in the
%           dummy variables of the kernels of f, and of their delta groups,
%           over the coefficients of all distributed monomials of f, see
%           kernel_stats. Fields, each a 1 x N array over the N spatial
%           variables f.vars.in unless stated otherwise:
%   dreq:   largest exponent that must be in the basis of Zop, over the
%           monomials with a numeric coefficient (nonzero constant term of
%           the kernel, not only decision variables): the largest second-
%           largest exponent of the free variables of a monomial, or the
%           largest free exponent if a delta group is present;
%   e1:     largest exponent over the monomials without delta group. At
%           most one free variable of a term of g has an exponent above d,
%           at most d+2p+1, as integration over the primary variable of s^k
%           adds the integral s^(k+1) of one limit of integration;
%   eM:     largest exponent of a delta group, s^k with k<=2p;
%   e1req, eMreq: as e1, eM, over the monomials with a numeric coefficient;
%   D:      largest total degree of a monomial;
%   hasD, hasNL: whether delta groups, and free variables, appear. By PSD
%           of the Gram matrix, the multiplier term of Zop is only needed
%           if hasD (delta part of g = square of terms with multiplier, so
%           zero if f has none), and the integral terms only if hasNL (if
%           all terms of f are delta terms, the multiplier-only Gram block
%           reproduces f);
%   ngrp:   largest number of delta groups in one term;
%   gmax:   largest number of dummy variables in one delta group (0 if
%           there is none);
%   hasdvar:  scalar, whether any kernel has decision variables;
%   hastie:   scalar, whether a term has two or more delta groups in some
%             variable and a group of three or more in some variable.
%   Terms with a zero kernel are skipped; terms whose kernel has only
%   decision variables count. A dummy variable evaluated at a boundary
%   point (a or b) is not in a group; it sets hasD.
%
% NOTES
% - otype='auto' gives the smallest type whose terms have the structure of
%   the terms of f. This is necessary, not sufficient: a certificate may 
%   need rows whose terms cancel in g. If F>=0 is infeasible for the otype 
%   chosen, try otype 3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_ineq_opts
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/22/2026: Initial coding, moved from 'piesos_ineq'; option
%                   otype='auto';

if ~isa(f,'polyopvar')
    error("Distributed polynomial functional must be specified as object of type 'polyopvar'.")
end
if ~f.chkval
    f = initialize(f);
end
if nargin<2 || isempty(opts)
    opts = struct();
end

% Statistics of the exponents and delta groups of the kernels of f, per
% spatial variable, over the coefficients of all distributed monomials (see
% OUTPUTS and kernel_stats), and hasdvar, whether f has decision variables
nvar = numel(f.vars.in);
z0 = zeros(1,nvar);
kst = struct('dreq',z0,'e1',z0,'e1req',z0,'eM',z0,'eMreq',z0,'D',z0, ...
               'hasD',false(1,nvar),'hasNL',false(1,nvar),'hasdvar',false, ...
                 'ngrp',z0,'gmax',z0,'hastie',false);
for i=1:size(f.degmat,1)
    Cop = f.C.ops{1,i};     % f.C is 1 x nZ, column i acts on monomial i
    if isa(Cop,'double')
        % For degree-0 monomials, degree of spatial variables is 0
        continue
    elseif isa(Cop,'dpvar') || isa(Cop,'polynomial')
        % For degree-0 monomials, no spatial variables should be present
        Cop = combine(Cop);
        if numel(Cop.varname)>0
            error("Enforcing positivity of polynomial vector fields is not supported.")
        end
    elseif isa(Cop,'tensopvar') && ~any(f.C.depmat.in(i,:))
        % Constants are stored as 'tensopvar' (see constant2polyopvar),
        % with no dependence on spatial variables, so no integration
    elseif isa(Cop,'intop')
        % Exponents of kernels and delta terms per variable (kernel_stats)
        kst = kernel_stats(kst,Cop,f.vars.in);
    else
        error("Coefficients acting on distributed monomials must be specified as 'intop' objects.")
    end
end
% Type of output of Zop for otype='auto', from the delta groups, set before
% the basis is built
if isfield(opts,'otype') && ~isnumeric(opts.otype) && strcmpi(opts.otype,'auto')
    opts.otype = ineq_otype(kst);
end
% Distributed monomial basis of g, for the resolved otype
Z = ineq_monomials(f,opts);
if isempty(Z)
    return
end
% Do not overwrite the degrees and exclusions specified by the user
[degs,exclude] = ineq_degs(f,Z,opts,kst);
% The degrees are set here, unless specified in opts.degs, or in opts.deg
if ~isfield(opts,'degs') && ~isfield(opts,'deg')                            % DJ, 11/07/2026
    opts.degs = degs;
end
if ~isfield(opts,'exclude') && ~isempty(exclude)
    opts.exclude = exclude;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function kst = kernel_stats(kst,Cop,varnames)
% KST = KERNEL_STATS(KST,COP,VARNAMES) updates the statistics KST of the
% exponents of the monomials in the dummy variables of the kernels of the
% scalar functional COP ('intop'), along each spatial variable (names
% VARNAMES in f). The dummy variables of a variable that are evaluated at
% the same point (intmat entries [p,p], a delta term) form one group, with
% the exponents added, and are the primary variable of Zop (multiplier
% term). The other dummy variables are free.
%
% INPUTS
% - kst:      struct of statistics, see the OUTPUTS of get_ineq_opts;
% - Cop:      'intop' object, the coefficient of one distributed monomial;
% - varnames: 1 x N cellstr, the spatial variables of f, to which the
%             fields of kst refer;
%
% OUTPUTS
% - kst:      the statistics, updated with the terms of Cop.

Kc = kernels(Cop);
Nc = numel(Cop.vars.in);
[~,gv] = ismember(Cop.vars.in,varnames);    % index of the variables of Cop in f
dnames = Cop.varnames.in_dum;               % names of the dummy variables
dv = Cop.depmat.in;                         % number of dummy variables
szK = size(Kc);
if Nc==1
    szK = [numel(Kc),1];
end
for L=1:numel(Kc)
    Kk = Kc{L};
    % Monomials of the kernel (zero kernels are skipped), and which have a
    % numeric coefficient: nonzero constant term (first row of C, for dpvar)
    if isa(Kk,'dpvar') || isa(Kk,'polynomial')
        if nnz(Kk.C)==0
            continue
        end
        degmat = full(Kk.degmat);       vn = Kk.varname;
        is_req = true(max(size(degmat,1),1),1);
        if isa(Kk,'dpvar')
            C1 = Kk.C(1,:);                 % constant term (first row)
            kst.hasdvar = kst.hasdvar || nnz(Kk.C)>nnz(C1);
            is_req = any(reshape(full(C1~=0),size(degmat,1),[]),2);
        end
    else
        if nnz(Kk)==0
            continue
        end
        degmat = zeros(1,0);            vn = cell(0,1);
        is_req = true;
    end
    if isempty(degmat)
        degmat = zeros(1,0);
    end
    % Term of Cop (row of intmat) for each spatial variable
    t = cell(1,max(Nc,1));
    [t{:}] = ind2sub(szK,L);
    % Number and largest size of the delta groups of this term
    ngrp_t = 0;     gmax_t = 0;
    for v=1:Nc
        if dv(v)==0
            continue
        end
        V = gv(v);
        % Dummy variable j is integrated over [s_A(j),s_B(j)], evaluated at
        % s_A(j) if A(j)==B(j): group with variable A(j). A(j) = 0 or dv+1
        % is a boundary point, which g cannot reach (a point evaluation)
        pr = Cop.intmat{v}(t{v},:);
        A = pr(1:2:end);    B = pr(2:2:end);
        grp = 1:dv(v);      is_bnd = false(1,dv(v));
        for j=find(A==B)
            if A(j)>=1 && A(j)<=dv(v) && A(j)~=j
                grp(grp==grp(j)) = grp(A(j));
            else
                is_bnd(j) = true;
            end
        end
        if any(is_bnd)
            kst.hasD(V) = true;     % as before: multiplier
        end
        % Exponents of the dummy variables in each monomial, added per group
        Ev = zeros(max(size(degmat,1),1),dv(v));
        for j=find(~is_bnd)
            c = find(strcmp(vn,dnames{v}{j}));
            if ~isempty(c) && ~isempty(degmat)
                Ev(:,j) = degmat(:,c);
            end
        end
        [ug,~,gi] = unique(grp(~is_bnd));
        if isempty(ug)
            continue
        end
        gsz = accumarray(gi(:),1)';
        jfree = find(~is_bnd);
        Ge = zeros(size(Ev,1),numel(ug));
        for k=1:numel(ug)
            Ge(:,k) = sum(Ev(:,jfree(gi==k)),2);
        end
        isD = gsz>=2;
        kst.hasD(V) = kst.hasD(V) || any(isD);
        kst.hasNL(V) = kst.hasNL(V) || any(~isD);
        % Number and size of the delta groups (points), for ineq_otype
        kst.ngrp(V) = max(kst.ngrp(V),nnz(isD));
        kst.gmax(V) = max([kst.gmax(V),gsz(isD)]);
        ngrp_t = max(ngrp_t,nnz(isD));
        gmax_t = max([gmax_t,gsz(isD)]);
        % Exponent that must be in the basis of Zop (second-largest, or largest
        % free exponent for a delta group), per monomial, and total degree
        kst.D(V) = max(kst.D(V),max(sum(Ge,2)));
        if any(isD)
            ed = zeros(size(Ge,1),1);
            if any(~isD)
                ed = max(Ge(:,~isD),[],2);
            end
            kst.eM(V) = max(kst.eM(V),max(max(Ge(:,isD))));
            if any(is_req)
                kst.eMreq(V) = max(kst.eMreq(V),max(max(Ge(is_req,isD))));
            end
        else
            Ge = sort(Ge,2,'descend');
            ed = zeros(size(Ge,1),1);
            if size(Ge,2)>=2
                ed = Ge(:,2);
            end
            kst.e1(V) = max(kst.e1(V),max(Ge(:,1)));
            if any(is_req)
                kst.e1req(V) = max(kst.e1req(V),max(Ge(is_req,1)));
            end
        end
        if any(is_req)
            kst.dreq(V) = max(kst.dreq(V),max(ed(is_req)));
        end
    end
    % Two points in a variable and three or more factors at one point: a
    % partial tie of copies (see ineq_otype)
    kst.hastie = kst.hastie || (ngrp_t>=2 && gmax_t>=3);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function otype = ineq_otype(kst)
% OTYPE = INEQ_OTYPE(KST) returns the type of output of Zop for
% otype='auto' from the delta groups of the terms of f, by the rule of the
% NOTES of get_ineq_opts, with a warning for a partial tie.
%
% INPUTS
% - kst:    struct, the statistics of the kernels of f, see the OUTPUTS of
%           get_ineq_opts (fields ngrp, gmax, hastie used);
%
% OUTPUTS
% - otype:  1, 2 or 3, the option otype of `piesos_quadopvar`.

% Tensor-product rows are needed for two points in a variable, pointwise
% rows for three or more factors at one point
is_tens = any(kst.ngrp>=2);
is_pointwise = any(kst.gmax>=3);
if ~is_tens
    otype = 1;
elseif ~is_pointwise
    otype = 2;
else
    otype = 3;
end
if kst.hastie
    warning('PIESOS:ineq:partialTie',['Some terms of f tie a subset of the copies of a',...
        ' variable, which no otype represents on the whole domain: f>=0 may be',...
        ' infeasible unless opts.intmat includes regions with these copies tied.']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Z = ineq_monomials(f,opts)
% Z = INEQ_MONOMIALS(F,OPTS) returns the distributed monomial basis of g:
% the monomials of half the degrees of those of F within half its Newton
% polytope, without those mixing states with different spatial dependence.
% With tensor-product rows, each monomial of different states also in
% each other order of its factors (dmonomials with 'full'): the pointwise
% rows of otype 3 couple only with the first factor of a tensor-product
% row, so that with one order the cone depends on the order of the states.
% piesos_quadvar keeps the pointwise rows of one order (opts.unique).
%
% INPUTS
% - f:      'polyopvar' object, the scalar functional of get_ineq_opts;
% - opts:   struct, the options of piesos_ineq with otype resolved
%           (numeric): tensor-product rows if otype, or that of opts.L or
%           opts.R, is 2 or 3 (or not given: the default 1 of
%           piesos_sosvar);
%
% OUTPUTS
% - Z:      nZ x 1 'polyopvar' object, the basis, with the sizes and
%           spatial dependence of the states of f; [] if no monomial
%           remains.

% Determine the degrees of the distributed monomials that appear in f
f_degmat = f.degmat;
if isempty(f_degmat)
    error("The distributed polynomial contains no monomials.")
end

% Initialize a set of monomial degrees to parameterize g
maxdeg = full(max(sum(f_degmat,2)));     % maximum total degree of the monomials
mindeg = full(min(sum(f_degmat,2)));     % minimum total degree of the monomials (for matrixvars, this will be at least 2)
Z_degmat = monomials(size(f_degmat,2),(floor(mindeg/2):ceil(maxdeg/2)));

% Discard unnecessary monomials
maxdegree = sparse(max(f_degmat,[],1)/2);      % row of max degrees in each variable
mindegree = sparse(min(f_degmat,[],1)/2);      % row of min degrees in each variable
Zdummy1 = bsxfun(@minus,maxdegree,Z_degmat);   % maxdegree monomial minus each monomial
Zdummy2 = bsxfun(@minus,Z_degmat,mindegree);   % each monomial minus mindegree monomial
[I,~] = find([Zdummy1 Zdummy2]<0);             % rows which contain negative terms
IND = setdiff(1:size(Z_degmat,1),I,'stable');  % rows not listed in I
Z_degmat = Z_degmat(IND,:);                    % discard all monomials rows listed in I
if ~isempty(Z_degmat)
    % Get rid of monomials outside of Newton polytope
    Z2 = f_degmat/2;
    Z_degmat = inconvhull(full(Z_degmat),full(Z2));
    Z_degmat = makesparse(Z_degmat);
    % Discard monomials of states with different spatial dependence: the
    % degree of a monomial in each variable is the sum over its factors of
    % the dependence of the states, and piesos_sosvar errors if these differ
    % along the variables on which it depends
    dep_mon = Z_degmat*double(f.depmat.in);         % nZ x N
    is_mixed = false(size(Z_degmat,1),1);
    for i=1:size(Z_degmat,1)
        is_mixed(i) = numel(unique(nonzeros(dep_mon(i,:))))>1;
    end
    % A certificate of f>=0 may need the square of a discarded monomial,
    % e.g. f = int (x1*x2)^2 for x1(s1), x2(s1,s2): f>=0 is then
    % conservative, and f=0 is imposed by piesos_ineq if none remain
    if any(is_mixed)
        msg = '%d of %d monomials in the basis for g mix states';
        msg = [msg,' depending on different variables (not supported)'];
        msg = [msg,' and are discarded: f>=0 may be conservative.'];
        if all(is_mixed)
            msg = [msg,' No monomials remain: f=0 is imposed.'];
        end
        warning('PIESOS:ineq:mixedDependence',msg, ...
            nnz(is_mixed),numel(is_mixed));
    end
    Z_degmat = Z_degmat(~is_mixed,:);
end
if isempty(Z_degmat)
    Z = [];
    return
end

% Declare the distributed monomial basis
% Pass dims and depmat, so states keep their size and spatial dependency
vartab = polyopvar(f.varname,f.dims.in,f.vars.in,f.dom.in,f.depmat.in);
[~,vidx] = ismember(vartab.varname,f.varname);
% Degrees per state: a state of size n is one variable of dmonomials. With
% tensor-product rows, each monomial in each order of its factors ('full')
otype_dflt = 1;     % the default of piesos_sosvar
if isfield(opts,'otype')
    otype_dflt = opts.otype;
end
is_tens = false;
for side = 'LR'
    otype_side = otype_dflt;
    if isfield(opts,side) && isfield(opts.(side),'otype')
        otype_side = opts.(side).otype;
    end
    is_tens = is_tens || any(otype_side(:,1)>=2);
end
if is_tens
    Z = dmonomials(vartab,Z_degmat(:,vidx),'full');
else
    Z = dmonomials(vartab,Z_degmat(:,vidx));
end
% dmonomials leaves monomials with zero coefficient in Z.Z, e.g. x1*x2 for
% the basis [x1^2;x2^2], which piesos_sosvar would treat as monomials:
% remove
Z = clean_poly(Z);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [degs,exclude] = ineq_degs(f,Z,opts,kst)
% [DEGS,EXCLUDE] = INEQ_DEGS(F,Z,OPTS,KST) returns the monomial degrees and
% term types of Zop for the basis Z, from the kernel statistics KST of F
% (see the NOTES of piesos_ineq).
%
% INPUTS
% - f:      'polyopvar' object, the scalar functional of get_ineq_opts;
% - Z:      'polyopvar' object, the basis of g;
% - opts:   struct, the options of piesos_ineq (field psatz used);
% - kst:    struct, the statistics of the kernels of f;
%
% OUTPUTS
% - degs:   struct with fields 'in' (nf x N) and 'out' (1 x N), see the
%           OUTPUTS of get_ineq_opts; 0 if Z has no spatial variables;
% - exclude: 1 x q cell, the option exclude of monopmat per block of Zop;
%           [] if Z has no spatial variables.

% Degrees and term types follow from the kernel exponents of f per spatial
% variable (see NOTES of piesos_ineq): d for the dummy variables of Zop
% (same for each factor), p for the primary variables, term types per
% block of Zop
[~,zvidx] = ismember(Z.C.vars.in,f.vars.in);
nfctr = max([1,max(Z.C.depmat.in(:))]);       % max number of factors
% Dummy degree d: the monomials with a numeric coefficient must be
% recovered (dreq). If f has decision variables, the other monomials may be
% set to zero by the solver, and a certificate may need others than those
% of f: then d, p are at least dbal, the smaller of
% ceil(D/4), D the maximal total degree of the kernels, and d=p such that
% e1<=d+2p+1 and eM<=2p
dbal = double(kst.hasdvar)*min(ceil(kst.D(zvidx)/4), ...
    max(ceil((kst.e1(zvidx)-1)/3),ceil(kst.eM(zvidx)/2)));
dz = max(kst.dreq(zvidx),dbal);
% Primary degree p: Pmat must have degree DP=max(eM,e1-d-1) in s, for the
% numeric monomials (req), which, for psatz multipliers of order q, is
% q+2*(p-ceil(q/2)) for p>=ceil(q/2)
DP = max([kst.eMreq(zvidx);kst.e1req(zvidx)-dz-1;zeros(1,numel(zvidx))],[],1);
qz = zeros(1,numel(zvidx));
if isfield(opts,'psatz')
    ps = opts.psatz;
    if islogical(ps)
        if isscalar(zvidx) && isvector(ps)
            ps = ps(:);
        end
        qz = size(ps,1:numel(zvidx))-1;
    elseif isscalar(ps)
        qz(:) = ps;
    elseif numel(ps)==numel(zvidx)
        qz = reshape(ps,1,[]);
    end
end
pz = ceil(DP/2);
pz(qz>0) = ceil(qz(qz>0)/2) + ceil(max(DP(qz>0)-qz(qz>0),0)/2);
pz = max(pz,dbal);
if isempty(zvidx)
    degs = 0;
    exclude = [];
else
    degs = struct('in',repmat(dz,nfctr,1),'out',pz);
    % Zop has a block per row of Z.C, acting on its output: types per row
    [~,zvidx_out] = ismember(Z.C.vars.out,f.vars.in);
    exclude = exclude_types(Z.C.depmat.out,[kst.hasD(zvidx_out);kst.hasNL(zvidx_out);kst.hasNL(zvidx_out)]);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function excl = exclude_types(depmat,allowed)
% EXCL = EXCLUDE_TYPES(DEPMAT,ALLOWED) returns the option 'exclude' of
% monopmat for each block (row of DEPMAT): a 3 x ... x 3 logical array over
% the variables on which the block depends, true for the combinations of
% types (multiplier, lower integral, upper integral) in which a type not
% allowed in that variable appears. A variable without allowed types (no
% terms in f) allows all, so that no block is empty.
%
% INPUTS
% - depmat:  q x N array, the dependence of the q blocks of Zop on the N
%            spatial variables (nonzero if a block depends on a variable);
% - allowed: 3 x N logical array, the term types allowed per variable;
%
% OUTPUTS
% - excl:    1 x q cell, element j the exclude array of block j (false for
%            a block without spatial variables).

q = size(depmat,1);
excl = cell(1,q);
for j=1:q
    idcs = find(depmat(j,:));
    if isempty(idcs)
        excl{j} = false;
        continue
    end
    A = 1;      nv = 0;
    for k=idcs
        al = double(allowed(:,k));
        if ~any(al)
            al = ones(3,1);
        end
        nv = nv+1;
        A = reshape(A(:)*al(:)',[3*ones(1,nv),1]);      % first variable fastest
    end
    excl{j} = ~logical(A);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

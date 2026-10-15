function [prog,g,Qmat,Zop,Z] = piesos_ineq(prog,f,opts)
% [PROG,G,QMAT,ZOP,Z] = PIESOS_INEQ(PROG,F,OPTS) adds the distributed SOS
% constraint F>=0, for a distributed polynomial F, to the PIESOS program
% PROG.
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program (see `piesos_program`);
% - f:      'polyopvar' object representing a distributed polynomial.
%           The object must be a scalar-valued functional: its
%           coefficients must be 'intop' objects (or constants for
%           degree-0 monomials), with f.matdim=[1,1] and no output
%           variables, as e.g. obtained using `innerprod`;
% - opts:   (optional) struct, passed on to `piesos_sosvar`, with fields
%               'psatz':  1 x N array of orders, passed on to
%                         `piesos_poslpivar` (see there);
%               'degs', 'exclude':  set here, based on the kernels of f
%                         (see NOTES), only if not specified;
%               'tol':    tolerance below which coefficients of f are
%                         taken to be zero. Defaults to 1e-14;
%
% OUTPUTS
% - prog:   'struct' specifying the same PIESOS program as the input, but
%           with the added constraint F>=0. This constraint is enforced
%           by declaring a distributed SOS polynomial G using
%           `piesos_sosvar`, and enforcing F==G using `piesos_eq`;
% - g:      'polyopvar' object representing the distributed SOS variable
%           for which F==G is enforced to impose positivity of F;
% - Qmat:   m x m cell of 'dpvar' blocks representing the symmetric
%           positive semidefinite (SOS in the primary variables)
%           multiplier operator (Pmat in `piesos_sosvar`) used to
%           parameterize g in the quadratic form
%               g = <Zop*Z(x), Qmat*Zop*Z(x)>;
% - Zop:    n x m 'tensopmat' object representing the monomial basis
%           operator used to parameterize g in the quadratic form;
% - Z:      m x 1 'polyopvar' object representing the distributed monomial
%           basis used to define g in the quadratic form;
%
% NOTES
% - If F is a 'dpvar', 'polynomial', or 'double' object, F>=0 is
%   enforced using `sosineq`, and if F is an 'opvar', 'opvar2d', 'dopvar',
%   or 'dopvar2d' object, using `lpi_ineq`. In these cases, only the
%   output PROG is returned.
%
% - A monomial in the basis for g (Z) cannot mix states that depend on
%   different spatial variables, e.g. x1(s1)*x2(s1,s2), see
%   `piesos_quadvar`. Such monomials are discarded with a warning (id
%   PIESOS:ineq:mixedDependence). This makes the constraint F>=0
%   conservative, in particular for F with terms such as (x1*x2)^2 that
%   require them. If no monomials remain, F=0 is imposed, with g=0 and
%   with Qmat, Zop and Z empty.
%
% - The degrees d (of the dummy variables in Zop) and p (of the primary
%   variables, in Qmat), and the types of terms in Zop, follow from the
%   exponents of the monomials in the dummy variables of the kernels of
%   f, per spatial variable. For a monomial, let e(1)>=e(2)>=... be the
%   sorted exponents of its free dummy variables, where dummy variables
%   evaluated at the same point (a delta term) count as one variable.
%   The terms of g, with dummy degree d and primary degree p, have all
%   but one exponent <=d, and that one <=d+2p+1 (from int_{s_i}^b s^k
%   ds). So, to recover a monomial,
%       d>=e(2),    d+2p+1>=e(1),    2p>=e_delta,
%   where e_delta is the exponent of a delta term.
%
% - The multiplier term in Zop (exclude(1)) is included only along
%   variables in which f has delta terms: otherwise its Gram blocks are
%   zero, as the delta part of g is a square. The integral terms
%   (exclude(2:3)) are included only along variables in which f has free
%   dummy variables: otherwise the multiplier block reproduces f.
%
% - If f has no decision variables, d and p are the smallest degrees
%   recovering all monomials. The total degree of the kernels does not
%   give d: for f=(int s^3*x ds)^2 it is 6, but d=3 is needed.
%
% - If f has decision variables, only the monomials with a numeric
%   coefficient need to be recovered: the others may be set to zero by
%   the solver. However, a certificate may need other monomials than
%   those of f: for -dV of a Lyapunov functional V, d=5 and p=2,
%   recovering all monomials (e(1)=9), is infeasible, whereas d=p=3 is
%   feasible, so that the smallest d and p depend on the data, not on
%   the monomials. Therefore, d and p are at least those for the numeric
%   monomials, and at least the smaller of ceil(D/4), for D the maximal
%   total degree of the kernels, and the d=p for which all monomials
%   fit: d+2p+1>=e(1), 2p>=e_delta.
%
% - A psatz multiplier of order q limits the degree of Qmat in s to
%   q+2(p-ceil(q/2)).
%
% - The same d is used for each factor of a distributed monomial, as the
%   order of the factors in the kernels of f is not unique. If F>=0 is
%   infeasible for the d and p chosen, increase them using opts.degs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_ineq
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
% DJ, 10/15/2026: Initial coding;

% % Check the inputs
% Check the program structure
if ~isa(prog,'struct')
    error("PIESOS program structure must be specified as 'struct'.")
end
% Check the polynomial functional
if isa(f,'dpvar') || isa(f,'polynomial') || isa(f,'double')
    % Enforce positivity of finite-dimensional polynomials using sos_ineq;
    if nargout>1
        error("Multiple output arguments is not supported in this case.")
    end
    prog = sosineq(prog,f);
    return
elseif isa(f,'dopvar') || isa(f,'dopvar2d') || isa(f,'opvar') || isa(f,'opvar2d')
    % Enforce positivity of operator variables using lpi_ineq
    if nargout>1
        error("Multiple output arguments is not supported in this case.")
    end
    prog = lpi_ineq(prog,f);
    return
elseif ~isa(f,'polyopvar')
    error("Distributed polynomial functional must be specified as object of type 'polyopvar'.")
end
% Make sure f is properly declared, and a scalar-valued functional, as
% piesos_eq(prog,f-g) compares it with the scalar functional g
if ~f.chkval
    f = initialize(f); 
end
if ~isequal(f.matdim,[1,1]) || ~isempty(f.vars.out)
    error("Distributed polynomial must be a scalar-valued functional (no output variables).")
end
% Check the options
if nargin==2
    opts = struct();
elseif nargin==3
    if ~isa(opts,'struct')
        error("Options for enforcing positivity must be specified as struct.")
    end
else
    error("Too many input arguments.")
end
% Tolerance below which coefficients are taken to be zero
if isfield(opts,'tol') && ~isempty(opts.tol)
    tol = opts.tol;
else
    tol = 1e-14;
end
% Merge duplicate rows, drop zero rows, and remove unused variables
f = clean_poly(f,tol);

% Determine the degrees of the distributed monomials that appear in f
f_degmat = f.degmat;
if isempty(f_degmat)
    error("The distributed polynomial contains no monomials.")
end
% Determine the degrees of the finite-dimensional monomials appearing in
% the operators acting on each distributed monomial
% Statistics of the exponents of the monomials in the dummy variables of
% the kernels, per spatial variable, see kernel_stats; and hasdvar, whether
% f has decision variables
nvar = numel(f.vars.in);
z0 = zeros(1,nvar);
kst = struct('dreq',z0,'e1',z0,'e1req',z0,'eM',z0,'eMreq',z0,'D',z0, ...
    'hasD',false(1,nvar),'hasNL',false(1,nvar),'hasdvar',false);
for i=1:size(f_degmat,1)
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
    % conservative, and f=0 is imposed below if none remain
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
    % The only way for the polynomial to be positive semidefinite, is for
    % it to be zero
    g = constant2polyopvar(0);
    prog = enforce_eq(prog,f,g,tol);
    Qmat = cell(0,0);
    Zop = tensopmat();
    Z = polyopvar();
    return
end

% Declare the distributed monomial basis
% Pass dims and depmat, so states keep their size and spatial dependency
vartab = polyopvar(f.varname,f.dims.in,f.vars.in,f.dom.in,f.depmat.in);
[~,vidx] = ismember(vartab.varname,f.varname);
Z = dmonomials(vartab,Z_degmat(:,vidx));
% dmonomials leaves monomials with zero coefficient in Z.Z, e.g. x1*x2 for  
% the basis [x1^2;x2^2], which piesos_sosvar would treat as monomials: 
% remove
Z = clean_poly(Z);

% Declare distributed SOS functional
% Do not overwrite the degrees and exclusions specified by the user
% Degrees and term types follow from the kernel exponents of f per spatial
% variable (see NOTES): d for the dummy variables of Zop (same for each
% factor), p for the primary variables, term types per block of Zop
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
        if numel(zvidx)==1 && isvector(ps),     ps = ps(:);     end
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
if ~isfield(opts,'degs')
    if isempty(zvidx)
        opts.degs = 0;
    else
        opts.degs = struct('in',repmat(dz,nfctr,1),'out',pz);
    end
end
if ~isfield(opts,'exclude') && ~isempty(zvidx)
    opts.exclude = exclude_types(Z.C.depmat.in,[kst.hasD(zvidx);kst.hasNL(zvidx);kst.hasNL(zvidx)]);
end
[prog,g,Qmat,Zop] = piesos_sosvar(prog,Z,opts);

% Enforce f==g>=0
prog = enforce_eq(prog,f,g,tol);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function prog = enforce_eq(prog,f,g,tol)
% PROG = ENFORCE_EQ(PROG,F,G,TOL) adds the constraint F=G to the program, for F
% and G scalar 'polyopvar' functionals. The constant (degree-0) terms are
% equated directly, as they are stored as 'tensopvar' in F (see
% constant2polyopvar) but as 'intop' in G (see piesos_sosvar), and
% 'tensopvar' + 'intop' is not supported.
% The remaining terms are equated using piesos_eq.

[f,cf] = split_constant(f);
[g,cg] = split_constant(g);
if ~isempty(cf) || ~isempty(cg)
    if isempty(cf),  cf = 0;  end
    if isempty(cg),  cg = 0;  end
    c = cf-cg;
    if isnumeric(c)
        c = polynomial(c);
    end
    prog = soseq(prog,c);
end

has_f = size(f.Z,1)>0;
has_g = size(g.Z,1)>0;
if has_f && has_g
    prog = piesos_eq(prog,f-g,tol);
elseif has_f
    prog = piesos_eq(prog,f,tol);
elseif has_g
    prog = piesos_eq(prog,-g,tol);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [f,c] = split_constant(f)
% [F,C] = SPLIT_CONSTANT(F) removes the degree-0 monomial from the scalar
% 'polyopvar' F, and returns its coefficient C (empty if F has none), a
% 'double' (for the 'tensopvar' of constant2polyopvar), or a 'dpvar'
% (for the 'intop' of piesos_sosvar).

c = [];
is_const = ~any(f.Z,2);
for i=reshape(find(is_const),1,[])
    Cop = f.C.ops{1,i};
    if isa(Cop,'intop') && isscalar(Cop.params) && isempty(Cop.intmat)
        Kc = kernels(Cop);
        ci = Kc{1};
    elseif isa(Cop,'tensopvar') && isequal(size(Cop.ops),[1,1]) && isa(Cop.ops{1,1},'sopvar') && isscalar(Cop.ops{1,1}.params)
        ci = full(Cop.ops{1,1}.params{1});
        if isempty(ci)
            ci = 0;
        end
    elseif isnumeric(Cop)
        ci = Cop;
    else
        error("Constant coefficients must be scalar 'tensopvar' or 'intop' objects without variables.")
    end
    if isempty(c)
        c = ci;
    else
        c = c+ci;
    end
end
% Remove the constant monomial and its coefficient
if any(is_const)
    keep = ~is_const;
    f.C.ops = f.C.ops(:,keep);
    f.C.dims.in = f.C.dims.in(keep);
    f.C.depmat.in = f.C.depmat.in(keep,:);
    f.Z = f.Z(keep,:);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function kst = kernel_stats(kst,Cop,varnames)
% KST = KERNEL_STATS(KST,COP,VARNAMES) updates the statistics KST of the
% exponents of the monomials in the dummy variables of the kernels of the
% scalar functional COP ('intop'), along each spatial variable (names VARNAMES
% in f). The dummy variables of a variable that are evaluated at the same
% point (intmat entries [p,p], a delta term) form one group, with the
% exponents added, and are the primary variable of Zop (multiplier term). The
% other dummy variables are free. Per spatial variable, fields of KST:
% - dreq:   largest exponent that must be in the basis of Zop, over the
%           monomials with a numeric coefficient (nonzero constant term of
%           the kernel, not only decision variables): the largest second-
%           largest exponent of the free variables of a monomial, or the
%           largest free exponent if a delta group is present;
% - e1:     largest exponent over the monomials without delta group. At most
%           one free variable of a term of g has an exponent above d, at most
%           d+2p+1, as integration over the primary variable of s^k adds
%           the integral s^(k+1) of one limit of integration;
% - eM:     largest exponent of a delta group, s^k with k<=2p;
% - e1req, eMreq: as e1, eM, over the monomials with a numeric coefficient;
% - D:      largest total degree of a monomial;
% - hasD, hasNL: whether delta groups, and free variables, appear. By PSD of
%           the Gram matrix, the multiplier term of Zop is only needed if
%           hasD (delta part of g = square of terms with multiplier, so zero
%           if f has none), and the integral terms only if hasNL (if all terms
%           of f are delta terms, the multiplier-only Gram block reproduces f);
% - hasdvar: (not per variable) whether any kernel has decision variables.

Kc = kernels(Cop);
Nc = numel(Cop.vars.in);
[~,gv] = ismember(Cop.vars.in,varnames);    % index of the variables of Cop in f
dnames = Cop.varnames.in_dum;               % names of the dummy variables
dv = Cop.depmat.in;                         % number of dummy variables
szK = size(Kc);
if Nc==1,   szK = [numel(Kc),1];    end
for L=1:numel(Kc)
    Kk = Kc{L};
    % Monomials of the kernel (zero kernels are skipped), and which have a
    % numeric coefficient: nonzero constant term (first row of C, for dpvar)
    if isa(Kk,'dpvar') || isa(Kk,'polynomial')
        if nnz(Kk.C)==0,    continue,   end
        degmat = full(Kk.degmat);       vn = Kk.varname;
        is_req = true(max(size(degmat,1),1),1);
        if isa(Kk,'dpvar')
            C1 = Kk.C(1,:);                 % constant term (first row)
            kst.hasdvar = kst.hasdvar || nnz(Kk.C)>nnz(C1);
            is_req = any(reshape(full(C1~=0),size(degmat,1),[]),2);
        end
    else
        if nnz(Kk)==0,      continue,   end
        degmat = zeros(1,0);            vn = cell(0,1);
        is_req = true;
    end
    if isempty(degmat),     degmat = zeros(1,0);    end
    % Term of Cop (row of intmat) for each spatial variable
    t = cell(1,max(Nc,1));
    [t{:}] = ind2sub(szK,L);
    for v=1:Nc
        if dv(v)==0,    continue,   end
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
        if any(is_bnd),     kst.hasD(V) = true;     end     % as before: multiplier
        % Exponents of the dummy variables in each monomial, added per group
        Ev = zeros(max(size(degmat,1),1),dv(v));
        for j=find(~is_bnd)
            c = find(strcmp(vn,dnames{v}{j}));
            if ~isempty(c) && ~isempty(degmat),     Ev(:,j) = degmat(:,c);  end
        end
        [ug,~,gi] = unique(grp(~is_bnd));
        if isempty(ug),     continue,   end
        gsz = accumarray(gi(:),1)';
        jfree = find(~is_bnd);
        Ge = zeros(size(Ev,1),numel(ug));
        for k=1:numel(ug)
            Ge(:,k) = sum(Ev(:,jfree(gi==k)),2);
        end
        isD = gsz>=2;
        kst.hasD(V) = kst.hasD(V) || any(isD);
        kst.hasNL(V) = kst.hasNL(V) || any(~isD);
        % Exponent that must be in the basis of Zop (second-largest, or largest
        % free exponent for a delta group), per monomial, and total degree
        kst.D(V) = max(kst.D(V),max(sum(Ge,2)));
        if any(isD)
            ed = zeros(size(Ge,1),1);
            if any(~isD),   ed = max(Ge(:,~isD),[],2);  end
            kst.eM(V) = max(kst.eM(V),max(max(Ge(:,isD))));
            if any(is_req), kst.eMreq(V) = max(kst.eMreq(V),max(max(Ge(is_req,isD))));  end
        else
            Ge = sort(Ge,2,'descend');
            ed = zeros(size(Ge,1),1);
            if size(Ge,2)>=2,   ed = Ge(:,2);   end
            kst.e1(V) = max(kst.e1(V),max(Ge(:,1)));
            if any(is_req), kst.e1req(V) = max(kst.e1req(V),max(Ge(is_req,1))); end
        end
        if any(is_req),     kst.dreq(V) = max(kst.dreq(V),max(ed(is_req)));    end
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function excl = exclude_types(depmat,allowed)
% EXCL = EXCLUDE_TYPES(DEPMAT,ALLOWED) returns the option 'exclude' of
% monopmat for each block (row of DEPMAT): a 3 x ... x 3 logical array over the
% variables on which the block depends, true for the combinations of types
% (multiplier, lower integral, upper integral) in which a type not allowed in
% that variable (3 x N logical array ALLOWED) appears. A variable without
% allowed types (no terms in f) allows all, so that no block is empty.

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
        if ~any(al),    al = ones(3,1); end
        nv = nv+1;
        A = reshape(A(:)*al(:)',[3*ones(1,nv),1]);      % first variable fastest
    end
    excl{j} = ~logical(A);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Test_polyopvar_subs_fctr.m -- substitution
% fnew = subs_factor(f,fctr,val), replacing factor number fctr by val,
% for a polynomial functional f (intop coefficient, no output variables).
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - Structural check only: fnew must be a validly declared 'polyopvar';
%   its value is not compared with an independent evaluation (see
%   Test_polyopvar_subs_const* for value checks).

% % Set parameters defining f and val
nterms = 1;     % number of terms in val
nfctrs = 1;     % number of factors in polynomial to substitute
d1 = 2;         % degree of monomial defining f
d2 = 1;         % maximal degree of polynomial val
nx1 = 1;        % number of state variables in f
nx2 = 1;        % number of state variables in val
nx_share = min(nx1,nx2);  % number of state varaibles in val that also appear in f
m = 1;          % maximal vector size of the states
Ntot = 1;       % number of spatial variables on which f depends

dnsty = 0.25;

% Set parameters defining the functional
otype = [2,2];
nr = 1;         % number of rows in f
use_full_int = false;
pdegv = 3;

% % Declare the state variable names
xidcs_full = randperm(nx1+nx2,nx1+nx2);
xidcs1 = sort(xidcs_full(1:nx1));
xidcs2 = sort(xidcs_full(nx1+1:end));
xidcs_full = [xidcs1; xidcs2];
xshare_idcs1 = sort(randperm(nx1,nx_share));
xshare_idcs2 = sort(randperm(nx2,nx_share));
xidcs2(xshare_idcs2) = xidcs1(xshare_idcs1);
xvarname_full = cell(nx1+nx2,1);
% All nx1+nx2 entries must be filled: xidcs1/xidcs2 index into them by values up to nx1+nx2.
for j=1:nx1+nx2
    xvarname_full{j} = ['x',num2str(xidcs_full(j))];
end
xvarname1 = xvarname_full(xidcs1);
xvarname2 = xvarname_full(xidcs2);

% % Set the size of the state variables
dims_full = randi(m,nx1+nx2,1);
dims1 = dims_full(xidcs1);
dims2 = dims_full(xidcs2);

% % Declare spatial variable names and domains
svarname_full = cell(1,Ntot);
for j=1:Ntot
    svarname_full{j} = ['s',num2str(j)];
end
dom_full = randi([0,1],[Ntot,1]) +[zeros(Ntot,1),rand([Ntot,1])];

% % Set dependency arrays for the variables
depmat_full = logical(randi([0,1],nx1+nx2,Ntot));
depmat1 = depmat_full(xidcs1,:);
depmat2 = depmat_full(xidcs2,:);

% % Set the variables and domains actually appearing in f and val
isvar1 = any(depmat1,1);                isvar2 = any(depmat2,1);
svars1 = svarname_full(isvar1);         svars2 = svarname_full(isvar2);
doms1 = dom_full(isvar1,:);             doms2 = dom_full(isvar2,:);
depmat1 = depmat1(:,isvar1);            depmat2 = depmat2(:,isvar2);
%degs1 = sum(depmat1,1);
pdegs = randi(2,[1,size(depmat2,2)]);
pdegsB = struct('in',randi(2,[1,size(depmat1,2)]),'out',randi(2,[1,size(depmat2,2)]));



% % Declare the polynomial function f
% Set the monomial defining f
fZ = randi([1,nx1],1,d1);
depmat_mid = depmat1(fZ,:);
% degs1 (degree of each spatial variable in the monomial Aop acts on, see
% rand_intop) must count repeated occurrences of a state in fZ; sum(depmat1,1)
% would build Aop for too low a degree and trip the input/output space check
% in @intop/mtimes_tensopvar.
degs1 = sum(depmat_mid,1);
dims_mid = reshape(dims1(fZ),[],1);
m1 = prod(dims1(fZ));
% Generate the functional acting on the monomial
Aop = rand_intop([nr,m1],svars1,doms1,degs1,pdegs);
% Generate f
f = polyopvar();
f.Z = fZ;
f.C = tensopmat(Aop);
f.varname = xvarname1;
f.dims.in = dims1;          f.dims.out = nr;
f.vars.in = svars1;         f.vars.out = cell(1,0);
f.dom.in = doms1;           f.dom.out = zeros(0,2);
f.depmat.in = depmat1;      f.depmat.out = zeros(1,0);



% % Declare the polynomial val
fctrs = randperm(d1,nfctrs);
val = [];
Bop_cell = cell(nfctrs,1);
vars = struct('in',{svars2},'out',{svars1});
dom = struct('in',doms2,'out',doms1);
depmat = struct('in',depmat2,'out',zeros(1,0));
for k=1:numel(fctrs)
    % Determine the function space of factor k
    idx = fZ(fctrs(k));
    depmat_k = depmat;
    depmat_k.out = depmat1(idx,:);
    m_k = dims1(idx);
    % Declare a new monomial to substitute
    Z_k = randi([1,nx2],[1,d2]);
    dims_k_2 = dims2(Z_k);
    depmat_k.in = depmat2(Z_k,:);
    % Declare a tensopvar object mapping monomials to factor
    dims_k_1 = round(m_k/d2)*ones(d2,1);
    if otype(1)==2 && prod(dims_k_1,1)~=m_k
        dims_k_1(1) = m_k;
        dims_k_1(2:end) = 1;
    elseif otype(1)==1
        dims_k_1 = m_k*ones(d2,1);
    end
    dims_k = [dims_k_1,dims_k_2];
    Bop_k = rand_tensopvar(dims_k,vars,dom,depmat_k,nterms,otype,pdegsB);
    Bop_cell{k} = Bop_k;
    % Declare the polynomial for substitution
    val_k = polyopvar();
    val_k.Z = Z_k;
    val_k.C = tensopmat(Bop_k);
    val_k.varname = xvarname2;
    val_k.dims.in = dims2;          val_k.dims.out = m_k;
    val_k.vars.in = svars2;         val_k.vars.out = svars1;
    val_k.dom.in = doms2;           val_k.dom.out = doms1;
    val_k.depmat.in = depmat1;      val_k.depmat.out = depmat_k.out;
    val_k.chkval = true;
    % Concatenate
    val = [val; val_k];
end


% % Perform the substitution
fnew = subs(f,fctrs,val);
if ~isa(fnew,'polyopvar')
    error("Test_polyopvar_subs_fctr failed: substitution did not return a 'polyopvar' object.")
elseif ~fnew.chkval
    error("Test_polyopvar_subs_fctr failed: substitution result is not marked as validly-declared.")
end
disp('Test_polyopvar_subs_fctr passed (structural check only -- see note above).')
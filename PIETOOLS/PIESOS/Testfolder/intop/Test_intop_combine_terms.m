% Test_intop_combine_terms.m -- Tests @intop/combine_terms.m: reordering
% the dummy variables of a functional acting on a distributed monomial of
% nx states over N=2 spatial variables must not change its action. Oracle:
% apply_op of the functional before vs. after combine_terms, on the same
% random monomial.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - nx=2 (degree-2 states, di=4 in each direction) is used because
%   rand_intop cannot generate the functional for nx=3 and depmat=true
%   (di=6 in each direction): split_ints then expands each term into up to
%   6!=720 sub-terms per direction, exceeding the maximal array size.
% - Value check only (relative residual); the canonical form of the result
%   is tested in Test_intop_combine_terms_canonical.m.

% Set parameters defining the distributed monomial
N = 2;                      % number of spatial variables
nx = 2;                      % number of state variables
dmin = 2;
dmax = 2;                   % maximal distributed monomial degree
degs = randi([dmin,dmax],1,nx);    % degree of each state in the monomial
d = sum(degs);
state_idcs = repelem(1:nx,degs);
n_arr = randi(2,[1,nx]);     % row dimension of state variables
q_arr = ones(1,nx);          % column dimensions of state variables (should be 1)
depmat = true(nx,N);
n_arr_in = n_arr(state_idcs);
q_arr_in = q_arr(state_idcs);
depmat_in = depmat(state_idcs,:);

% Set parameters defining functional
m = 2;  
n = prod(n_arr_in);            % column dimension of functional       
dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
deps = sum(depmat_in,1);
pdegs = 3;

% Set parameters defining the function
q = 2;              % column dimension of function
pdegv = 2;

% Generate a random intop object
Aop = rand_intop([m,n],N,dom,deps,pdegs);
Aop2 = Aop;

% Reorder the variables
Anew = combine_terms(Aop2,degs,n_arr,depmat);

% Generate a random distributed monomial on which the functional acts
vnames = cell(1,N);
for i=1:N
    vnames{i} = ['s',num2str(i)];
end
xcell = cell(1,nx);
pvars_arr = cell(1,nx);
for i=1:nx
    dep_i = depmat(i,:);
    vnames_i = vnames(dep_i);
    pvars_arr{i} = polynomial(vnames_i(:))';
    xcell{i} = rand_poly([n_arr(i),q_arr(i)],pvars_arr{i},pdegv,50);
    if all(isequal(xcell{i},0))
        xcell{i} = eye(n_arr(i),q_arr(i));
    end
end
vnames_in = Aop.varnames.in;
ycell = cell(1,d);
xval = 1;
cntr = zeros(1,N);
pvars_new = cell(1,d);
for j=1:d
    dep_j = depmat_in(j,:);
    vnames_j = vnames_in(dep_j);
    cntr(dep_j) = cntr(dep_j)+1;    % move to the next variable
    vars_j = cellfun(@(a,b) a(b),vnames_j,num2cell(cntr(dep_j)),'UniformOutput',false);
    pvars_new{j} = polynomial([vars_j{:}])';
    ycell{j} = subs(xcell{state_idcs(j)},pvars_arr{state_idcs(j)},pvars_new{j});
    xval = kron(xval,ycell{j});
end

% Apply the functionals
Ax = apply_op(Aop2,xval);
Ax_alt = apply_op(Anew,xval);


% Compare
tol = 1e-8;
err_A = cleanpoly(Ax - Ax_alt,1e-10);
if ~isdouble(err_A)
    error("At least one of the outputs is not constant.")
elseif norm(double(err_A),'fro')/norm(double(Ax),'fro')>tol
    error("Outputs do not match: relative residual exceeds tolerance.")
else
    disp('err_C = ')
    disp(double(err_A)/norm(double(Ax)))
end
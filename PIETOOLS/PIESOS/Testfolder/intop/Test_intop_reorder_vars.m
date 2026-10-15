% Test_intop_reorder_vars.m -- Tests reorder_vars for 'intop' (N=1).
%
% Oracle: for a random functional Aop on one spatial variable with d dummy
% copies, Anew = reorder_vars(Aop,state_order) applied to the monomial with
% its factors in the order state_order (dummy variables renamed to the slot
% each factor now occupies) must give the same value as Aop applied to the
% original monomial, both evaluated with apply_op.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Only N=1 is tested (depmat is all ones, one dummy copy per factor).

% Set parameters defining the distributed monomial
N = 1;                      % number of variables (TEST ONLY MAKES SENSE FOR N=1)
d = 4;                      % distributed monomial degree
n_arr = randi(1,[1,d]);     % row dimension of state variables
p_arr = ones(1,d);          % column dimensions of state variables (should be 1)
depmat = logical(ones(d,N));    % which state variable depends on which spatial var

% Set parameters defining functional
m = 1;  
n = prod(n_arr);            % column dimension of functional       
dom = [0,1];
degs = sum(depmat,1);
pdegs = 3;

state_order = randperm(degs,degs);  % Order of the variables

% Set parameters defining the function
p = 2;              % column dimension of function
pdegv = 2;

% Generate random intop objects
Aop = rand_intop([m,n],N,dom,degs,pdegs);

% Generate a random distributed monomial on which the functional acts
vnames = Aop.varnames.in;
xcell = cell(1,d);
xval = 1;
cntr = zeros(1,d);
pvars_arr = cell(1,d);
for i=1:d
    dep_i = depmat(i,:);
    vnames_i = vnames(dep_i);
    cntr(dep_i) = cntr(dep_i)+1;    % move to the next variable
    vars_i = cellfun(@(a,b) a(b),vnames_i,num2cell(cntr(dep_i)),'UniformOutput',false);
    pvars_arr{i} = polynomial([vars_i{:}])';
    xcell{i} = rand_poly([n_arr(i),p_arr(i)],pvars_arr{i},pdegv,50);
    if all(isequal(xcell{i},0))
        xcell{i} = eye(n_arr(i),p_arr(i));
    end
    xval = kron(xval,xcell{i});
end

% Reorder the variables
Anew = reorder_vars(Aop,state_order);

% Reorder the distributed monomial accordingly
depmat_new = depmat(state_order,:);
vnames_new = Aop.varnames.in;
xcell_new = xcell(state_order);
xval_new = 1;
cntr = zeros(1,d);
pvars_new = cell(1,d);
for i=1:d
    dep_i = depmat_new(i,:);
    vnames_i = vnames_new(dep_i);
    cntr(dep_i) = cntr(dep_i)+1;    % move to the next variable
    vars_i = cellfun(@(a,b) a(b),vnames_i,num2cell(cntr(dep_i)),'UniformOutput',false);
    pvars_new{i} = polynomial([vars_i{:}])';
    xcell_new{i} = subs(xcell_new{i},pvars_arr{state_order(i)},pvars_new{i});
    xval_new = kron(xval_new,xcell_new{i});
end

% Apply the functionals
Av = apply_op(Aop,xval);
Av_alt = apply_op(Anew,xval_new);


% Compare: residual must be constant and, relative to |Av_alt|, below tol
tol = 1e-8;
err_A = cleanpoly(Av - Av_alt,1e-10);
nrmref = norm(full(Av_alt.C),'fro');
if ~isdouble(err_A)
    error("At least one of the outputs is not constant.")
elseif norm(double(err_A),'fro')>tol*max(1,nrmref)
    error("Outputs do not match: residual exceeds tolerance.")
else
    disp('err_C = ')
    disp(double(err_A))
end

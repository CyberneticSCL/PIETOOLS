function fval = apply_functional(Kcell,xvals,degmat,idx_mat,vars,dom)
% FVAL = APPLY_FUNCTIONAL(KCELL,XVALS,IDX_MAT,VARS,DOM) computes the value
% of the integral
% f(x) = sum_{i=1}^{m} int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
%           KCELL{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1
% where kj = idx_ii(i,j) for j=1,...,d and i=1,...,m for m, for the
% given values of the polynomial functions xi;
%
% INPUTS
% - Kcell:  m x 1 cell with the ith element a polynomial of d variables,
%           specifying the kernel to be used in term i of the integral 
%           (functional) operator;
% - xvals:  n x 1 array of type 'polynomial' 
% - degmat: 1 x n array of integers, specifying the degrees of the
%           state variables appearing in the distributed monomials;
% - idx_mat:    m x d array specifying for each of the kernels the
%               associated order of the variables of integration, so that
%               if idx_mat(l,:) = [i,j,k], then a <= ti <= tj <= tk <= b
%               and therfore the lth term has integral
%                   int_{a}^{b} int_{ti}^{b} int_{tj}^{b} ... dtk dtj dti
%               where d=sum(degmat)
% - vars:   d x 1 array of type 'polynomial' specifying the dummy variables
%           for integration, (t1,...,td). These must be the same as the
%           variables that appear in the polynomials in Kcell!
% - dom:    1 x 2 array, [a,b], specifying the spatial domain
%
% OUTPUTS
% - fval:   The value of the function f(x) for the specified K and x;


d = sum(degmat);
nstates = size(degmat,2);
if numel(xvals)~=nstates
    error("Distributed monomial must be specified as 1 x d array")
end
if size(degmat,1)~=1
    error("Polynomials involving multiple monomials are not supported.")
end
% Establish for each factor in the monomial which state variable is
% considered
state_idcs = [];
for ii=1:nstates
    state_idcs = [state_idcs; ii*ones(degmat(ii),1)];
end
xvals_full = polynomial(zeros(d,1));
for ii=1:d
    state_num = state_idcs(ii);
    var_ii = xvals(state_num).varname;
    if isscalar(var_ii)
        xvals_full(ii) = xvals(state_num);
        xvals_full(ii).varname = vars(ii).varname;
    elseif isempty(var_ii)
        xvals_full(ii) = xvals(state_num);
    elseif numel(var_ii)>1
        error("Each state variable can depend on at most one independent variable")
    end
end

if isempty(idx_mat)
    % List all possible orders of the variables t1 through td
    %   idx_mat(l,:) = [i,j,k] means a <= ti <= tj <= tk <= b in term l
    idx_mat = 1;
    for ii=2:d
        n_ords = size(idx_mat,1);
        idx_mat_new = zeros(ii*n_ords,ii);
        for jj=1:ii
            % Place variable t_ii in position jj
            idx_mat_new(jj:ii:end,:) = [idx_mat(:,1:jj-1),ii*ones(n_ords,1),idx_mat(:,jj:end)];
        end
        idx_mat = idx_mat_new;
    end
end


fval = zeros(size(Kcell{1}));
if d>2 && numel(Kcell)>factorial(d)
    error("Number of operators should be d! for monomial degree d.")
end
for ii=1:min(numel(Kcell),factorial(d))
    if isempty(Kcell{ii})
        continue
    end
    % Check the order of the variables:
    %   idx_ii = [i,j,k] implies a <= ti <= tj <= tk <= b
    idx_ii = idx_mat(ii,:);
    
    % Perform the appropriate integrals
    %   int_{a}^{b} int_{t_k1}^{b} ... int_{t_kd}^{b} 
    %           Kcell{i}(t1,...,td)*x1(t1)...xd(td) dt_kd ... dt_k1 
    % where kj = idx_ii(i,j) for j=1,...,d and i=1,...,m for m=d!.
    fval_ii = Kcell{ii};
    for jj=d:-1:2
        var_num = idx_ii(jj);   % integrating over the jth biggest variable
        L = vars(idx_ii(jj-1)); % integrating from the j-1th biggest variable up to b
        fval_ii = fval_ii*xvals_full(var_num);
        fval_ii = int(fval_ii,vars(var_num),L,dom(2));
    end
    % Finally, perform the integral from a to b in the smallest variable
    var_num = idx_ii(1); 
    fval_ii = fval_ii*xvals_full(var_num);
    fval_ii = int(fval_ii,vars(var_num),dom(1),dom(2));
    fval = fval + fval_ii;
end

if d==2 && numel(Kcell)==3
    % Third term corresponds to the multiplier
    %   int_{a}^{b} P(s)*x(s)^2 ds
    fval = fval + int(Kcell{3}*xvals_full(1)*subs(xvals_full(2),vars(2),vars(1)),vars(1),dom(1),dom(2));
end


end
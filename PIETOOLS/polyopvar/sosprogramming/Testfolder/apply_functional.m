function fval = apply_functional(Kcell,xvals,idx_mat,vars,dom)
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
% - xvals:  d x 1 array of type 'polynomial' 
% - idx_mat:    m x d array specifying for each of the kernels the
%               associated order of the variables of integration, so that
%               if idx_mat(l,:) = [i,j,k], then a <= ti <= tj <= tk <= b
%               and therfore the lth term has integral
%                   int_{a}^{b} int_{ti}^{b} int_{tj}^{b} ... dtk dtj dti
% - vars:   d x 1 array of type 'polynomial' specifying the dummy variables
%           for integration, (t1,...,td). These must be the same as the
%           variables that appear in the polynomials in Kcell!
% - dom:    1 x 2 array, [a,b], specifying the spatial domain
%
% OUTPUTS
% - fval:   The value of the function f(x) for the specified K and x;


d = size(idx_mat,2);
if numel(xvals)~=d
    error("Distributed monomial must be specified as 1 x d array")
end
for ii=1:d
    var_ii = xvals(ii).varname;
    if isscalar(var_ii)
        xvals(ii).varname = vars(ii).varname;
    elseif numel(var_ii)>1
        error("Each state variable can depend on at most one independent variable")
    end
end

fval = zeros(size(Kcell{1}));
for ii=1:numel(Kcell)
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
        fval_ii = fval_ii*xvals(var_num);
        fval_ii = int(fval_ii,vars(var_num),L,dom(2));
    end
    % Finally, perform the integral from a to b in the smallest variable
    var_num = idx_ii(1); 
    fval_ii = fval_ii*xvals(var_num);
    fval_ii = int(fval_ii,vars(var_num),dom(1),dom(2));
    fval = fval + double(fval_ii);
end


end
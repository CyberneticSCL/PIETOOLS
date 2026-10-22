function Kop = rand_intop(dims,vars,dom,degs,pdegs,q)
% KOP = RAND_INTOP(DIMS,VARS,DOM,DEGS,PDEGS,Q) generates a random
% 'intop' object representing a PI functional.
%
% INPUTS
% - dims:   1 x 2 array specifying the matrix dimensions [m,n] of the
%           functional. A scalar specifies m=n;
% - vars:   1 x N cellstr specifying the names of the spatial variables
%           along the N spatial directions.
%           ALTERNATIVELY: a scalar, 'N', specifying the number of spatial
%           variables, which are named {'s1',...,'sN'};
% - dom:    N x 2 array of integers, specifying the interval on which each
%           of the spatial variables is defined. A 1 x 2 array specifies
%           the same interval for each variable;
% - degs:   1 x N array specifying for each of the spatial variables the
%           degree of that variable in the distributed monomial on which
%           the functional acts. A scalar specifies the same degree for
%           each variable;
% - pdegs:  1 x N array specifying the degrees of the spatial variables in
%           the parameters (kernels) defining the functional. A scalar
%           specifies the same degree for each variable;
% - q:      (optional) nonnegative integer specifying the number of
%           decision variables 'coeff1',...,'coeffq' that the kernels
%           defining the functional depend on affinely. Defaults to 0, in
%           which case kernels are plain 'polynomial' objects. If q>0,
%           kernels are 'dpvar' objects (see rand_dpvar.m), all sharing the
%           same q decision variables;
%
% OUTPUTS
% - Kop:    m x n 'intop' object representing the randomly generated
%           functional;

% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/22/2026: Split the rows of intmat that are not ordered with
%   intop.canonicalize_intmat before declaring the functional, as the
%   'intop' constructor rejects them. split_ints follows the signed
%   convention for reversed limits, canonicalize_intmat the indicator;

% Check that the dimensions are appropriately specified
if isscalar(dims)
    dims = [dims,dims];
elseif numel(dims)>2
    error("Matrix-valued 'intop' object can have at most a row and column dimension.")
end
dims = reshape(dims,[1,2]);

% Check that the variables are appropriately specified
if isa(vars,'double') && isscalar(vars)
    % Declare N variables
    N = vars;
    vars = cell(1,N);
    for i=1:N
        vars{i} = ['s',num2str(i)];
    end
else
    % Otherwise, use variable names specified by the user
    N = numel(vars);
    if isa(vars,'polynomial')
        if ~ispvar(vars)
            error("Variables must be specified as 1xN array of pvar objects.")
        end
        vars = pvar2varname(vars(:))';
    elseif ~iscellstr(vars)
        error("Variables must be specified as 1xN 'cellstr' object.")
    end
end

% Check that the domain is properly specified
if size(dom,1)==1
    dom = repmat(dom,[N,1]);
elseif size(dom,1)~=N
    error("Number of spatial domains should match number of spatial variables.")
end

% Check that the degrees are appropriately specified
if isscalar(degs)
    degs = degs*ones(1,N);
elseif numel(degs)~=N
    error("Number of degrees should match number of variables.")
end
if any(degs<0) || any(round(abs(degs))~=degs)
    error("Degrees should be nonnegative.")
end

% Check that the degrees are appropriately specified
if isscalar(pdegs)
    pdegs = pdegs*ones(1,N);
elseif numel(pdegs)~=N
    error("Number of parameter degrees should match number of variables.")
end
if any(pdegs<0) || any(round(abs(pdegs))~=pdegs)
    error("Degrees should be nonnegative.")
end

% Check that the number of decision variables is appropriately specified
if nargin<6 || isempty(q)
    q = 0;
elseif ~isscalar(q) || q<0 || round(q)~=q
    error("Number of decision variables must be a nonnegative integer.")
end

% Declare spatial variables defining the kernels
pvars = polynomial(zeros(0,1));
for i=1:N
    % Declare the input variables
    di = degs(i);
    pvars_i = cell(di,1);
    for j=1:di
        pvars_i{j} = [vars{i},'_dum',num2str(j)];
    end
    pvars = [pvars; polynomial(pvars_i)];
end

% Declare a random intmat structure
intmat = cell(1,N);
pdegs_cell = cell(1,N);
for i=1:N
    di = degs(i);
    ntrms_i = min(randi([factorial(di+1),factorial(2*di)]),25);
    intmat_i = zeros(ntrms_i,2*di);
    intmat_i(:,1:2:end) = randi([0,di],[ntrms_i,di]);
    intmat_i(:,2:2:end) = randi([1,di+1],[ntrms_i,di]);
    intmat_i = uniquerows_integerTable(intmat_i);
    intmat{i} = intmat_i;
    pdegs_cell{i} = pdegs(i)*ones(1,di);
end
% Get rid of inconsistent integrals
intmat = remove_badints(intmat);
for i=1:N
    if size(intmat{i},1)==0
        intmat{i} = repmat([0,degs(i)+1],[1,degs(i)]);
    end
end
ntrms = cellfun(@(a)size(a,1),intmat);
pdegs = cell2mat(pdegs_cell)';

% Declare a random kernel for each integral
params = cell([ntrms,1]);
nZ_max = 50;               % maximal number of monomials in each parameter
if q==0
    for k=1:numel(params)
        params{k} = rand_poly(dims,pvars,pdegs,nZ_max);
    end
else
    % dpvar kernels: one shared dvarname list, needed for split_ints/intop
    % to combine the params consistently
    dvarname = cell(q,1);
    for i=1:q
        dvarname{i} = ['coeff',num2str(i)];
    end
    for k=1:numel(params)
        params{k} = rand_dpvar(dims,dvarname,pvars,pdegs,nZ_max);
    end
end

% Declare the functional
[intmat_new,params_new] = split_ints(intmat,params);
% Split the rows that are not ordered, repeating their kernels
try                                                                         % DJ (CC), 10/22/2026
    [intmat_new,params_new] = intop.canonicalize_intmat(intmat_new,params_new); % DJ (CC), 10/22/2026
catch                                                                       % DJ (CC), 10/22/2026
    % All regions are empty: use the full integral of the first kernel      % DJ (CC), 10/22/2026
    for i=1:N,  intmat_new{i} = repmat([0,degs(i)+1],1,degs(i));  end       % DJ (CC), 10/22/2026
    params_new = params_new(1);                                             % DJ (CC), 10/22/2026
end                                                                         % DJ (CC), 10/22/2026
Kop = intop(params_new,intmat_new,vars,dom);

end
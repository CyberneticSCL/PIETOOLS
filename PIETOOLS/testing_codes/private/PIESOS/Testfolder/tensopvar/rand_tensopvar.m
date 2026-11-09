function Pop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty)
% POP = RAND_TENSOPVAR(DIMS,VARS,DOM,DEPMAT,NTERMS,TYPE,DEGS,DNSTY)
% returns a 'tensopvar' object of specified dimensions and in the
% specified variables, with the factors in the tensor-PI operator declared
% as randomly generated 'sopvar' objects (rand_sopvar).
%
% INPUTS
% - dims:   d x 2 array specifying for each of the factors in the tensor-PI
%           operator the matrix dimensions of this factor;
% - vars:   struct with fields 'in' and 'out'. Each field is a cellstr
%           object, specifying the names of the input and output variables
%           of the tensor-PI operator. Can also be set to "vars.in=N" and
%           "vars.out=M" to randomly generate N input and M output
%           variables;
% - dom:    struct with fields 'in' and 'out'. Each field is a N x 2 (or
%           M x 2) array specifying for each of the variables on what
%           interval it is defined. A 1 x 2 array can be used for a field,
%           or for dom itself, to use the same interval for all variables;
% - depmat: struct with fields 'in' and 'out'. Each field is a d x N (or
%           d x M) logical array, specifying for each factor whether or not
%           that factor maps from (to) the specified input (output)
%           variables;
% - nterms: scalar integer specifying the number of terms in the tensor-PI
%           operator;
% - type:   (optional) 1 x 2 array with elements in {1,2}, set as the
%           'otype' of Pop: the type of product of the outputs and of the
%           inputs of the factors (1: pointwise product, 2: tensor
%           product). Defaults to [1,2];
% - degs:   (optional) struct with fields 'in' and 'out'. Each is a 1 x N
%           (or 1 x M) array of integers, specifying the maximal degree of
%           the monomials in the input and output variables in each of the
%           factors in the tensor-PI operators. Can also be set to a
%           scalar, for the same degree in all variables. Defaults to 3;
% - dnsty:  (optional) scalar in (0,1), specifying the density of the
%           sparse coefficient matrices used to generate the factors in
%           the tensor-PI operator. Defaults to [], for which rand_sopvar
%           chooses the density;
%
% OUTPUTS
% - Pop:    randomly generated 'tensopvar' object, with nterms terms of d
%           factors;

% DJ (CC), 10/15/2026: Initial coding;

% Set default values
if nargin<=5
    type = [1,2];
end
if nargin<=6
    degs = 3;
end
if nargin<=7
    dnsty = [];
end

% Check that the matrix dimensions are appropriate
d = size(dims,1);
if type(1)==1 && max(dims(:,1))~=min(dims(:,1))
    error("Row dimensions of operators should match for pointwise product of the outputs.")
end
if type(2)==1 && max(dims(:,2))~=min(dims(:,2))
    error("Column dimensions of operators should match for pointwise product of the inputs.")
end

% Check that the variables are properly specified
if ~isa(vars,'struct') || ~isfield(vars,'in') || ~isfield(vars,'out')
    error("Variables should be specified as struct with fields 'in' and 'out'.")
end
if isa(vars.in,'double')
    % Declare M=vars.in input variables
    N = vars.in;
    vars.in = cell(1,N);
    idcs = sort(randperm(2*N,N));
    for i=1:N
        vars.in{i} = ['s',num2str(idcs(i))];
    end
elseif ~iscellstr(vars.in)
    error("Input variables should be specified as 1 x N 'cellstr' object.")
end
vars.in = vars.in(:)';      N = numel(vars.in);
if isnumeric(vars.out) && isscalar(vars.out)
    % Declare M=vars.out output variables
    M = vars.out;
    vars.out = cell(1,M);
    idcs = sort(randperm(2*N+2*M,M));
    for i=1:M
        vars.out{i} = ['s',num2str(idcs(i))];
    end
elseif ~iscellstr(vars.out)
    error("Output variables should be specified as 1 x M 'cellstr' object.")
end
vars.out = vars.out(:)';    M = numel(vars.out);

% Check that the domains are appropriately specified
if isnumeric(dom) && all(size(dom)==[1,2])
    dom1 = dom;
    dom = struct();
    dom.in = repmat(dom1,N,1);
    dom.out = repmat(dom1,M,1);
elseif ~isa(dom,'struct') || ~isfield(dom,'in') || ~isfield(dom,'out')
    error("Domains should be specified as struct with fields 'in' and 'out'.")
end
if all(size(dom.in)==[1,2])
    % Assume same domain for all input variables
    dom.in = repmat(dom.in,N,1);
elseif ~all(size(dom.in)==[N,2])
    error("Input domains should be specified as N x 2 array for N input variables.")
end
if all(size(dom.out)==[1,2])
    % Assume same domain for all output variables
    dom.out = repmat(dom.out,M,1);
elseif ~all(size(dom.out)==[M,2])
    error("Output domains should be specified as M x 2 array for M output variables.")
end

% Check that the dependency arrays are appropriate
if ~isa(depmat,'struct') || ~isfield(depmat,'in') || ~isfield(depmat,'out')
    error("'depmat' should be specified as struct with fields 'in' and 'out'.")
elseif size(depmat.in,1)~=d || size(depmat.out,1)~=d
    error("Row dimensions of 'depmat' should match number of factors.")
elseif size(depmat.in,2)~=N || size(depmat.out,2)~=M
    error("Column dimensions of 'depmat' should match number of variables.")
elseif ~islogical(depmat.in) || ~islogical(depmat.out)
    error("'depmat' arrays should be of type 'logical'");
end

% Check that the degrees are properly specified
if isnumeric(degs) && isscalar(degs)
    deg1 = degs;
    degs = struct();
    degs.in = repmat(deg1,1,N);
    degs.out = repmat(deg1,1,M);
elseif ~isa(degs,'struct') || ~isfield(degs,'in') || ~isfield(degs,'out')
    error("Monomial degrees should be specified as struct with fields 'in' and 'out'.")
end
if isscalar(degs.in)
    % Assume same degree for all input variables
    degs.in = repmat(degs.in,1,N);
elseif numel(degs.in)~=N
    error("Number of input monomial degrees should match number of input variables.")
end
if isscalar(degs.out)
    % Assume same degree for all output variables
    degs.out = repmat(degs.out,1,M);
elseif numel(degs.out)~=M
    error("Number of output monomial degrees should match number of output variables.")
end

% Declare a tensor-PI operator with randomly generated factors
Pop = tensopvar();
Pop.otype = type;
Pop.ops = cell(nterms,d);
for j=1:d
    matdim_j = dims(j,:);
    vars_j = vars;
    vars_j.in = vars.in(depmat.in(j,:));
    vars_j.out = vars.out(depmat.out(j,:));
    dom_j = dom;
    dom_j.in = dom.in(depmat.in(j,:),:);
    dom_j.out = dom.out(depmat.out(j,:),:);
    degs_j = degs;
    degs_j.in = degs.in(depmat.in(j,:));
    degs_j.out = degs.out(depmat.out(j,:));
    for i=1:nterms
        Pop.ops{i,j} = rand_sopvar(matdim_j,vars_j,dom_j,degs_j,dnsty);
    end
end
Pop = initialize(Pop);

end
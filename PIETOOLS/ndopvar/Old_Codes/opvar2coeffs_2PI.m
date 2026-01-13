function [C,dim,d] = opvar2coeffs(Pop,d)
% [C,D] = OPVAR2COEFSS(POP) takes an opvar object 'Pop' and returns a cell
% 'C' and scalar 'd' representing the parameters defining 'Pop' as
%   Pop.R.R1(s,th) = kron(eye(m),Z_{d}(s)) C{1} kron(eye(n),Z_{d}(th))
%   Pop.R.R2(s,th) = kron(eye(m),Z_{d}(s)) C{2} kron(eye(n),Z_{d}(th))
% where Z_{d}(s) is the basis of monomials of degree at most 'd' in s, and
% where dim=[m,n].
%
% INPUTS
% - Pop:    mxn 'opvar' object representing a 2-PI operator
% - d:      (optional) scalar 'double' specifiying a desired degree of the
%           monomial basis. Defaults to the maximal degree of the monomials
%           appearing in Pop.R.R1 and Pop.R.R2.
% 
% OUTPUTS
% - C:      1x2 'cell' specifying the coefficients representing Pop.R.R1
%           (C{1}) and Pop.R.R2 (C{2}) in the quadratic form.
% - dim:    1x2 'double' array specifying the dimensions [m,n] of the
%           opvar object.
% - d:      degree of the monomial basis used in the quadratic
%           representation.

if ~isa(Pop,'opvar')
    error("Input  must be of type 'opvar'.")
end
if any(Pop.dim(1,:))
    error("Only operators from L2 to L2 are supported.")
end
if ~isempty(Pop.R.R0) && any(any(~isequal(polynomial(Pop.R.R0),0)))
    error("Only 2-PI operators are supported.")
end
if nargin==1
    d = inf;
end

% Extract the relevant information
m = Pop.dim(2,1);       n = Pop.dim(2,2);
var1_name = Pop.var1.varname{1};
var2_name = Pop.var2.varname{1};
R1 = polynomial(Pop.R.R1);
R2 = polynomial(Pop.R.R2);

% Determine the maximal monomial degree in both var1 and var2
dmax = max([max(R1.degmat,[],"all"),max(R2.degmat,[],"all")]);
if nargin==1
    d = dmax;
elseif d<dmax
    error("Specified degree is smaller than maximal monomial degree of parameters.")
end

% Get coefficients
C1 = get_quadratic_form(R1,d,var1_name,var2_name);
C2 = get_quadratic_form(R2,d,var1_name,var2_name);
C = {C1,C2};

end



%%
function C = get_quadratic_form(R,d,var1_name,var2_name)
% Compute coefficients C such that R(x,y)=(Im o Zd(x))^T C (In o Zd(y)),
% where x is var1 and y is var2, and where o denotes the Kronecker product.
% INPUTS
% - R:      mxn 'polynomial' class object in variables var1 and/or var2, of
%           degree at most d in either variable (not cumulative)
% - d:      scalar 'double' specifiying the degree of the monomial basis
%           Zd in var1 and var2
% - var1:   char object specifying the name of the variable "x"
% - var2:   char object specifying the name of the variable "y"
%
% OUTPUTS
% - C:      m(d+1) x n(d+1) array of type 'double' specifying the
%           coefficients such that R(x,y)=(Im o Zd(x))^T C (In o Zd(y))

% Extract degrees and coefficients
[m,n] = size(R);
Rvars = R.varname;
degs = R.degmat;
coeffs = R.coeff;

if isempty(Rvars) || (isscalar(degs) && degs==0)
    % If R does not involve any variables, R=C;
    r_idcs = repmat(1:m,[1,n]);
    c_idcs = repmat(1:n,[m,1]);
    C = sparse(r_idcs(:),c_idcs(:),coeffs(:),m*(d+1),n*(d+1));
    return
elseif any(~ismember(Rvars,{var1_name;var2_name}))
    error("Parameters in the opvar can only depend on var1 and var2.")
end
if ~ismember(var1_name,Rvars)
    % Degrees are only in var2
    % --> add zero degrees for var1
    degs = [zeros(size(degs)),degs];
elseif ~ismember(var2_name,Rvars)
    % Degrees are only in var1
    % --> add zero degrees for var2
    degs = [degs, zeros(size(degs))];
elseif ~strcmp(var1_name,Rvars{1})
    % Make sure var1 is in the first column of degs
    degs = fliplr(degs);
end

% Set the coefficients
vals = [];
c_idcs = [];
r_idcs = [];
rplus = repmat((0:m-1)'.*(d+1),[1,n]);      rplus = rplus(:)';
cplus = repmat((0:n-1).*(d+1),[m,1]);       cplus = cplus(:)';
for deg2=0:d
    % Check which monomials are of degree deg2 in var2
    z_nums = degs(:,2)==deg2;
    nz = nnz(z_nums);
    if nz==0
        continue
    end
    % Extract nonzero coefficients
    vals_tmp = coeffs(z_nums,:);
    is_nonzero = vals_tmp(:)~=0;
    vals_tmp = reshape(vals_tmp(is_nonzero),[],1);
    vals = [vals; vals_tmp(is_nonzero)];
    % Establish associated rows in C
    r_idcs_tmp = repmat(degs(z_nums,1)+1,[1,m*n]);
    r_idcs_tmp = reshape(r_idcs_tmp + repmat(rplus,[nz,1]),[],1);
    r_idcs = [r_idcs; r_idcs_tmp(is_nonzero)];
    % Establish associated columns in C
    c_idcs_tmp = (deg2+1)*ones(nnz(z_nums),m*n);
    c_idcs_tmp = reshape(c_idcs_tmp + repmat(cplus,[nz,1]),[],1);
    c_idcs = [c_idcs; c_idcs_tmp(is_nonzero)];
end
C = sparse(r_idcs,c_idcs,vals,m*(d+1),n*(d+1));

end
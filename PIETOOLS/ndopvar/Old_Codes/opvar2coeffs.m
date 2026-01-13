function [C_str] = opvar2coeffs(Pop,d)
% [C_STR] = OPVAR2COEFSS(POP) takes an opvar object 'Pop' and returns a
% struct 'C_str' with fields C_str.C representing the parameters defining
%   Pop.R.R0(s)    = kron(eye(m),Z_{d}(s)) C{0}
%   Pop.R.R1(s,th) = kron(eye(m),Z_{d}(s)) C{1} kron(eye(n),Z_{d}(th))
%   Pop.R.R2(s,th) = kron(eye(m),Z_{d}(s)) C{2} kron(eye(n),Z_{d}(th))
% where Z_{d}(s) is the basis of monomials of degree at most d in s, for
% d = C_str.deg and [m,n]=C_str.dim;
%
% INPUTS
% - Pop:    mxn 'opvar' object representing a 2-PI operator
% - d:      (optional) scalar 'double' specifiying a desired degree of the
%           monomial basis. Defaults to the maximal degree of the monomials
%           appearing in Pop.R.R1 and Pop.R.R2.
% 
% OUTPUTS
% - C_str:  'struct' with fields:
%           C_str.deg - the degree 'd' of the basis of monomials used to
%           define the parameters Pop.R;
%           C_str.C - a 1x3 cell containing the coefficients representing
%           the parameters Pop.R in the quadratic format;
%           C_str.dim - the size of the operator Pop, Pop.dim(2,:);
%           C_str.dom - the spatial domain of the operator, Pop.I;
%           C_str.vars - 1x2 array of 'pvar' objects, representing s, theta

if ~isa(Pop,'opvar')
    error("Input  must be of type 'opvar'.")
end
if any(Pop.dim(1,:))
    error("Only operators from L2 to L2 are supported.")
end
if nargin==1
    d = inf;
end

% Extract the relevant information
dim = Pop.dim(2,:);
var1_name = Pop.var1.varname{1};
var2_name = Pop.var2.varname{1};
R0 = polynomial(Pop.R.R0);
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
C0 = get_quadratic_form(R0,d,var1_name,[]);
C1 = get_quadratic_form(R1,d,var1_name,var2_name);
C2 = get_quadratic_form(R2,d,var1_name,var2_name);
C = {C0,C1,C2};

% Collect parameters representing the operator in a struct
C_str = struct();
C_str.C = C;
C_str.dim = dim;
C_str.dom = Pop.I;
C_str.deg = d;
C_str.vars = [Pop.var1,Pop.var2];

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
    if isempty(var1_name) && isempty(var2_name)
        C = sparse(double(R));
    else
        r_idcs = repmat(1:m,[1,n]);
        c_idcs = repmat(1:n,[m,1]);
        if isempty(var2_name)
            % Express only in terms of var1
            C = sparse(r_idcs(:),c_idcs(:),coeffs(:),m*(d+1),n);
        elseif isempty(var1_name)
            C = sparse(r_idcs(:),c_idcs(:),coeffs(:),m,n*(d+1));
        else
            C = sparse(r_idcs(:),c_idcs(:),coeffs(:),m*(d+1),n*(d+1));
        end
    end
    return
%elseif any(~ismember(Rvars,{var1_name;var2_name}))
%    error("Parameters in the opvar can only depend on var1 and var2.")
end
if isempty(var2_name)
    % Express only in terms of var1
    rplus = repmat((0:m-1).*(d+1),[1,n]);
    r_idcs = degs(:,1)+1 + rplus;
    c_idcs = repmat(kron((1:n),ones(1,m)),[size(degs,1),1]);
    C = sparse(r_idcs(:),c_idcs(:),coeffs(:),m*(d+1),n);
    return
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
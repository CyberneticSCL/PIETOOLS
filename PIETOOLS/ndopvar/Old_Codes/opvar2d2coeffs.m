function [C_str] = opvar2d2coeffs(Pop,d)
% [C_STR] = OPVAR2COEFSS2D(POP) takes an opvar2d object 'Pop' and returns a
% struct 'C_str' with fields C_str.C representing the parameters defining
%   Pop.R22{1,1}(s)    = kron(eye(m),Z_{d}(s)) C{0}
%   Pop.R22{i,j}(s,th) = kron(eye(m),Z_{d}(s)) C{1} kron(eye(n),Z_{d}(th))
% where [m,n]=C_str.dim, d = C_str.deg, and 
%   Z_{d}(s) = kron(Z_{d}(s1),Z_{d}(s2))
% for Z_{d}(si) = [1;si;si^2;...;si^d] is the basis of monomials of degree 
% at most d in si.
%
% INPUTS
% - Pop:    mxn 'opvar2d' object representing a 9-PI operator
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

if ~isa(Pop,'opvar2d')
    error("Input  must be of type 'opvar2d'.")
end
Pop = poly_opvar2d(Pop);
if any(Pop.dim(1,:))
    error("Only operators from L2 to L2 are supported.")
end
if nargin==1
    d = inf;
end

% Extract the relevant information
dim = Pop.dim(4,:);
var1 = Pop.var1;
var2 = Pop.var2;

% Determine the maximal monomial degree in both var1 and var2
dmax = 0;
for ii=1:numel(Pop.R22)
    dmax = max(dmax,max(max(Pop.R22{ii}.degmat)));
end
if nargin==1
    d = dmax;
elseif d<dmax
    error("Specified degree is smaller than maximal monomial degree of parameters.")
end

% Get coefficients
C_cell = cell(size(Pop.R22));
C_cell{1} = get_quadratic_form(Pop.R22{1,1},var1,[],d);
C_cell{2,1} = get_quadratic_form(Pop.R22{2,1},var1,var2(1),d);
C_cell{3,1} = get_quadratic_form(Pop.R22{3,1},var1,var2(1),d);
C_cell{1,2} = get_quadratic_form(Pop.R22{1,2},var1,var2(2),d);
C_cell{1,3} = get_quadratic_form(Pop.R22{1,3},var1,var2(2),d);
for ii=[5,6,8,9]
    C_cell{ii} = get_quadratic_form(Pop.R22{ii},var1,var2,d);
end

% Collect parameters representing the operator in a struct
C_str = struct();
C_str.C = C_cell;
C_str.dim = dim;
C_str.dom = Pop.I;
C_str.deg = d;
C_str.vars = [Pop.var1,Pop.var2];

end
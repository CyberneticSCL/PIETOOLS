function [C,d] = opvar2coeffs(Pop,d)
% [C,D] = OPVAR2COEFSS(POP) takes an opvar object 'Pop' and returns a cell
% 'C' and scalar 'd' representing the parameters defining 'Pop' as
%   Pop.R.R1(s,th) = Z_{d}(s) C{1} Z_{d}(th)
%   Pop.R.R2(s,th) = Z_{d}(s) C{2} Z_{d}(th)
% where Z_{d}(s) is the basis of monomials of degree at most 'd' in s.
%
% INPUTS
% - Pop:    1x1 'opvar' object representing a 2-PI operator
% - d:      (optional) scalar 'double' specifiying a desired degree of the
%           monomial basis. Defaults to the maximal degree of the monomials
%           appearing in Pop.R.R1 and Pop.R.R2.
% 
% OUTPUTS
% - C:      1x2 'cell' specifying the coefficients representing Pop.R.R1
%           (C{1}) and Pop.R.R2 (C{2}) in the quadratic form.
% - d:      degree of the monomial basis used in the quadratic
%           representation

if ~isa(Pop,'opvar')
    error("Input  must be of type 'opvar'.")
end
if ~all(all(Pop.dim==[0,0;1,1]))
    error("Only scalar-valued operators from L2 to L2 are supported.")
end
if ~isempty(Pop.R.R0) && any(any(~isequal(Pop.R.R0,0)))
    error("Only 2-PI operators are supported.")
end
if nargin==1
    d = inf;
end

% Extract the relevant information
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
% Compute coefficients C such that R(x,y)=Zd(x)^T C Zd(y),
% where x is var1 and y is var2
% INPUTS
% - R:      1x1 'polynomial' class object in variables var1 and/or var2, of
%           degree at most d in either variable (not cumulative)
% - d:      scalar 'double' specifiying the degree of the monomial basis
%           Zd in var1 and var2
% - var1:   char object specifying the name of the variable "x"
% - var2:   char object specifying the name of the variable "y"
%
% OUTPUTS
% - C:      d+1 x d+1 array of type 'double' specifying the coefficients
%           such that R(x,y)=Zd(x)^T C Zd(y)

% Extract degrees and coefficients
Rvars = R.varname;
degs = R.degmat;
coeffs = R.coeff;

if isempty(Rvars)
    % If R does not involve any variables, R=C;
    C = sparse(d+1,d+1);
    C(1,1) = double(R);
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
for deg2=0:d
    % Check degrees of var1 for monomials of degree deg2 in var2
    r_nums = degs(:,2)==deg2;
    if ~any(r_nums)
        continue
    end
    r_idcs = [r_idcs; degs(r_nums,1)+1];
    c_idcs = [c_idcs; (deg2+1)*ones(nnz(r_nums),1)];
    vals = [vals; coeffs(r_nums)];
end
C = sparse(r_idcs,c_idcs,vals,d+1,d+1);

end
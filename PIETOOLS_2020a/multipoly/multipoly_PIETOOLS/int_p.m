function b = int_p(a,x,x1,x2)
%
% DESCRIPTION
%   A blatent copy of diff.m to
%   Integrate a polynomial.
%
% INPUTS
%   A: polynomial
%   X: integrate with respect to X.
%   x1: lower limit of integration
%   x2: upper limit of integration
%
% OUTPUTS
%   B: polynomial
%
% SYNTAX
%   B = int(A,X);
%     Integrate the polynomial, A, with respect to X, from x1 to x2.  A
%     should be a polynomial and X should be a polynomial variable or a
%     string. x1 and x2 should be scalars or scalar polynomials.
%     Integration is done element-by-element if A is a matrix.

% 6/20/2010: MMP  Initial Coding   -   based on diff.m by PJS
% 6/16/2013: MMP  Modified to eliminate for-loops

% Error Checking
if nargin~=4
    error('Error in calling diff');
end

if isa(x,'polynomial')
    x = x.varname{1};
elseif ~ischar(x)
    error('X must be a polynomial variable or a string');
end

% Get polynomial info about a
a = polynomial(a);
%a = combine(a);
adeg = a.degmat;
avar = a.varname;
nta = size(adeg,1);
[nra,nca]=size(a);
acoef = reshape(a.coefficient,nta,nra*nca);

% Get format x1 and x2
if ~isa(x1,'double');
    x1 = polynomial(x1);
    x1 = combine(x1);
end
if ~isa(x2,'double');
    x2 = polynomial(x2);
    x2 = combine(x2);
end


% Find variable we are integrating with respect to.
varnumb=find(strcmp(x,avar));  % MMP - 6.16.2013

if isempty(varnumb)
    b=a*(x2-x1);
else
    
    % Integrate
    %   for i1 = 1:nta
    %     acoef(i1,:) = acoef(i1,:)/(adeg(i1,varnumb)+1);
    %   end
    %     acoef = bsxfun(@times,acoef,1./(adeg(:,varnumb)+1));
    acoef = diag(1./(adeg(:,varnumb)+1))*acoef;
    adeg(:,varnumb) = adeg(:,varnumb)+1;
    a_int = polynomial(acoef,adeg,avar,[nra nca]);
    b = combine(subs_p(a_int,x,x2)-subs_p(a_int,x,x1));
end




function [prog,g,Qmat,Zop,Z] = piesos_ineq(prog,f,opts)
% [PROG,G,QMAT,ZOP,Z] = PIESOS_INEQ(PROG,F,OPTS) takes a PIESOS program
% PROG with distributed polynomial variable F and adds the distributed SOS 
% constraint F>=0 to the program.
%
% INPUTS
% - prog:   PIESOS program structure;
% - f:      'polyopvar' object representing a distributed polynomial
%           variable;
% - opts:   (optional) struct with field 'psatz', to be passed to
%           piesos_sosvar to indicate whether or not to include a psatz
%           term;
%
% OUTPUTS
% - prog:   lpiprogram structure specifying the same optimization program
%           as the input, but with the added constraint f >= 0. This
%           constraint is enforced by declaring a distributed SOS
%           polynomial, g, using piesos_sosvar, and enforcing f==g using
%           piesos_eq;
% - g:      'polyopvar' object representing the distributed SOS variable
%           for which f==g is enforced to impose positivity of f;
% - Qmat:   n x n 'dpvar' object representing the symmetric positive 
%           semidefinite multiplier operator used to parameterize g in the 
%           quadratic form,
%               g = <Zop*Z(x), Qmat*Zop*Z(x)>
% - Zop:    n x m 'tensopvar' object representing the monomial basis
%           operator used to parameterize g in quadratic form;
% - Z:      m x 1 'polyopvar' object representing the distributed monomial
%           basis used to define g in quadratic form;
%
% NOTES:
% The degrees of the finite-dimensional monomials used to parameterize the
% multiplier Qmat and operator Zop are determined based on the degrees of
% the finite-dimensional monomials encountered in the coefficient operators
% defining f. Here, a balance is aimed to be struck between numerical
% complexity and conservatism of the constraint, by trying to match the
% monomials that appear in f in a minimal manner. However, experiments
% suggest that even if all monomials that appear in the coefficients of f
% also appear in the coefficients of g, increasing this number of monomials
% may still reduce conservatism, at the cost of higher complexity.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_ineq
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 03/21/2026: Initial coding

% % Check the inputs
% Check the program structure
if ~isa(prog,'struct')
    error("PIESOS program structure must be specified as 'struct'.")
end
% Check the polynomial functional
if isa(f,'dpvar') || isa(f,'polynomial') || isa(f,'double')
    % Enforce positivity of finite-dimensional polynomials using sos_ineq;
    if nargout>1
        error("Multiple output arguments is not supported in this case.")
    end
    prog = sosineq(prog,f);
    return
elseif isa(f,'dopvar') || isa(f,'dopvar2d') || isa(f,'opvar') || isa(f,'opvar2d')
    % Enforce positivity of operator variables using lpi_ineq
    if nargout>1
        error("Multiple output arguments is not supported in this case.")
    end
    prog = lpi_ineq(prog,f);
    return
elseif ~isa(f,'polyopvar')
    error("Distributed polynomial functional must be specified as object of type 'polyopvar'.")
end
if numel(f.pvarname)>1
    error("Multivariate case is currently not supported.")
end
% Check the options
if nargin==2
    opts = struct();
elseif nargin==3
    if ~isa(opts,'struct')
        error("Options for enforcing positivity must be specified as struct.")
    end
else
    error("Too many input arguments.")
end

% Determine the degrees of the distributed monomials that appear in f
f_degmat = f.degmat;
if isempty(f_degmat)
    error("The distributed polynomial contains no monomials.")
end
% Determine the degrees of the finite-dimensional monomials appearing in
% the operators acting on each distributed monomial
excludeL = [1,0,0]';
pdegs_full = cell(size(f_degmat,1),1);
pdegs = 0;
for i=1:size(f_degmat,1)
    Cop = f.C.ops{i};
    if isa(Cop,'double')
        % For degree-0 monomials, degree of spatial variables is 0
        pdegs_full{i} = [0,0];
        continue
    elseif isa(Cop,'dpvar') || isa(Cop,'polynomial')
        % For degree-0 monomials, no spatial variables should be present
        Cop = combine(Cop);
        if numel(Cop.varname)>0
            error("Enforcing positivity of polynomial vector fields is not supported.")
        end
        pdegs_full{i} = [0,0];
    elseif isa(Cop,'intop')
        % Check if the functional includes a delta term
        if any(Cop.omat==0)
            excludeL(1) = 0;
        end
        % Check the maximal (cumulative) degree of the monomials
        C_degmat = Cop.params.degmat;
        if isempty(C_degmat)
            pdegs_full{i} = [0,0];
        else
            pdegs_full{i} = [max(max(C_degmat)),max(sum(C_degmat,2))];
        end
        % d = sum(f_degmat(i,:));
        % pdeg2 = ceil(pdegs_full{i}(2)/(4*d));       % degree of dummy variables
        % pdeg1 = ceil((pdegs_full{i}(1)-pdeg2)/2);   % degree of primary variable
        % pdegs = max(pdegs,[pdeg1,pdeg2]);
        pdegs = max(pdegs,ceil(pdegs_full{i}(2)/4));
    else
        error("Coefficients acting on distributed must be specified as 'intop' objects.")
    end
end

% Initialize a set of monomial degrees to parameterize g
maxdeg = full(max(sum(f_degmat,2)));     % maximum total degree of the monomials
mindeg = full(min(sum(f_degmat,2)));     % minimum total degree of the monomials (for matrixvars, this will be at least 2)
Z_degmat = monomials(size(f_degmat,2),(floor(mindeg/2):ceil(maxdeg/2))); 

% Discard unnecessary monomials
maxdegree = sparse(max(f_degmat,[],1)/2);      % row of max degrees in each variable
mindegree = sparse(min(f_degmat,[],1)/2);      % row of min degrees in each variable
Zdummy1 = bsxfun(@minus,maxdegree,Z_degmat);   % maxdegree monomial minus each monomial
Zdummy2 = bsxfun(@minus,Z_degmat,mindegree);   % each monomial minus mindegree monomial
[I,~] = find([Zdummy1 Zdummy2]<0);             % rows which contain negative terms
IND = setdiff(1:size(Z_degmat,1),I,'stable');  % rows not listed in I
if isempty(IND)
    % The only way for the polynomial to be positive semidefinite, is for
    % it to be zero
    prog = piesos_eq(prog,f-g);
    Qmat = 0;
    Zop = 0;
    Z = 0;
    return
else
    Z_degmat = Z_degmat(IND,:);                       % discard all monomials rows listed in I
end

% Get rid of monomials outside of Newton polytope
Z2 = f_degmat/2;
Z_degmat = inconvhull(full(Z_degmat),full(Z2));
Z_degmat = makesparse(Z_degmat);

% Declare the distributed monomial basis
vartab = polyopvar(f.varname,f.pvarname,f.dom); % fine in 1D, doesn't work if not all states depend on all variables
% vartab = f;
% vartab.degmat = zeros(0,numel(f.varname));
% vartab.C = tensopvar();
Z = dmonomials(vartab,Z_degmat);

% Declare distributed SOS functional
opts.deg = pdegs;
opts.exclude = excludeL;
[prog,g,Qmat,Zop] = piesos_sosvar(prog,Z,opts);

% Enforce f==g>=0
prog = piesos_eq(prog,f-g);

end
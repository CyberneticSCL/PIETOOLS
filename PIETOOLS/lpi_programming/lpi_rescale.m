function [prog_new,vartab_old2new] = lpi_rescale(prog,dom_new)
% PROG_NEW = LPI_RESCALE(PROG,DOM_NEW) takes an LPI program structure 
% 'prog' defined on some domain prog.dom, and performs a transformation
% of the independent variables to exist on the desired domain 'dom_new',
% returning a structure 'prog' that represents the same program but now
% with prog_new.dom = dom_new.
%
% INPUT
% - prog:       'struct' specifying an LPI optimization program structure
%               (see also 'lpiprogram'). The format should match that used
%               by SOSTOOLS;
% - dom_new:    nx2 array specifying the desired domain on which the
%               spatial variables in the optimization program should exist.
%               Defaults to [-1,1];
%
% OUTPUT
% - prog_new:   'struct' specifying the same optimization program as passed
%               by the user, but now with the spatial variables and 
%               monomial vectors updated to exist on the new domain, so
%               that prog_new.dom = dom_new;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_rescale
%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 10/23/2024

% % Check that the program structure is properly specified.
if ~isa(prog,'struct')
    error('LPI program structure should be specified as object of type ''struct''.')
elseif ~isfield(prog,'dom')
    error('The program structure does not include a spatial domain -- please use ''lpiprogram'' to initialize your program');
end
% Extract the domain of the spatial variables in the old structure.
dom_old = prog.dom;
nvars = size(prog.dom,1);

% % Check that the spatial domain is properly specified.
if nargin<2 || isempty(dom_new)
    % Default to [-1,1];
    dom_new = [-ones(nvars,1),ones(nvars,1)];
elseif ~isa(dom_new,'double') && ~(isa(dom_new,'polynomial') && isdouble(dom_new))
    error("Spatial domain should be specified as nx2 array of type 'double'.")
elseif size(dom_new,2)~=2
    error("Spatial domain should be specified as nx2 array.")
elseif size(dom_new,1)~=nvars
    if size(dom_new,1)==1
        % Assume same domain for all spatial variables.
        dom_new = repmat(dom_new,[nvars,1]);
    else
        error("Number of rows in domain array should match number of spatial variables in the original program structure.")
    end
end
dom_new = double(dom_new);
% Make sure the domain makes sense.
if any(dom_new(:,2)<=dom_new(:,1))
    error('Upper boundary of domain should be strictly greater than lower boundary.')
end

% Deal with case that the variables already exist on the desired domain.
if all(all(dom_old==dom_new))
    prog_new = prog;
    vars_new = prog.var;
    return
end


% % Define the variable substitution.

% First, extract the primary spatial variables, dummy variables, and free
% variables, the last of which do not need to exist on any domain.
vartab = prog.vartable;
var1 = vartab(1:nvars,1);         % primary variables
var2 = vartab(nvars+1:2*nvars,1); % dummy variables
varf = vartab(2*nvars+1:end,1);   % free variables, not constrained to spatial domain

% For each variable xi in [dom_old(i,1),dom_new(i,2)], introduce a new 
% variable yi in [dom_new(i,1),dom_new(i,2)] such that xi=dom_old(i,1) 
% implies yi=dom_new(i,1), and xi=dom_old(i,2) implies yi=dom_new(i,2).
a_old = dom_old(:,1);       L_old = dom_old(:,2)-dom_old(:,1);
a_new = dom_new(:,1);       L_new = dom_new(:,2)-dom_new(:,1);

% y(x) for x in [a_old, a_old+L_old]
var1_old2new = (var1-a_old).*(L_new./L_old) +a_new;
var2_old2new = (var2-a_old).*(L_new./L_old) +a_new;
vartab_old2new = [var1_old2new; var2_old2new; varf];

% x(y) for y in [a_new, a_new+L_new]
var1_new2old = (var1-a_new).*(L_old./L_new) +a_old;
var2_new2old = (var2-a_new).*(L_old./L_new) +a_old;
vartab_new2old = [var1_new2old; var2_new2old; varf];


% Substitute
prog_new = prog;
for ii=1:numel(prog_new.expr.Z)
    % Take the degrees 'degmat_full' of monomials for constraint ii.
    degmat_full = prog.expr.Z{ii};
    if ~any(any(degmat_full))
        % The monomials don't depend on the spatial variables
        % --> no need to change anything.
        continue
    end
    % Determine unique degrees 'degmat' and permutation matrix Pmat s.t.
    %   degmat_full = Pmat*degmat;
    [Pmat,degmat] = uniquerows_integerTable(degmat_full);
    % Build a vector of monomials for the given degrees;
    Zii_old = polynomial(prod((vartab').^degmat,2));
    % Express in terms of new variables (on new domain);
    Zii_new = polynomial(subs(Zii_old,vartab,vartab_new2old));
    % Determine permutation matrix Cii s.t.
    %   Zii_new = Cii*(vartab_new.^degmat_new)
    degmat_new = Zii_new.degmat; %NO, Zii.varname may no longer include all variables!
    Cii = Zii_new.C';
    % if size(C_ii,1)==size(C_ii,2) && ~any(any(C_ii-speye(size(C_ii))))
    %     % The new monomials are the same as the old ones
    %     % --> no need to change anything.
    %     continue
    % end

    % But how to express constraints in terms of new variables?
    % We have constraints
    %   (A*xi-b)*bkldiag(Z1,Z2,...)
    % for monomial vectors Z1, Z2, ...
    % But, SOSTOOLS stores only Z:=[Z1;Z2;...].
    % If we know that Z1; Z2 all start with the constant monomial 1, we can
    % determine where each new monomial vector Zi starts. But, it is not
    % strictly necessary for each Zi to include the constant monomial 1, so
    % we do not know where in Z any monomial vector Zi starts and ends...
    % Without that information, we cannot express Zi in terms of the new
    % monomials.
    
end


end


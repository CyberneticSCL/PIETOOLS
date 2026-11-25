function F = constant2polyopvar(C,vars,dom)
% F = CONSTANT2POLYOPVAR(C,VARS,DOM) takes a constant, a decision variable
% or a known function C of spatial variables, and converts it to a
% degree-0 distributed monomial F.
%
% INPUTS
% - C:      'double', 'polynomial', or 'dpvar' column vector (or
%           'polyopvar', returned unchanged). A 'polynomial' may depend on
%           spatial variables in VARS; a 'dpvar' cannot depend on any;
% - vars:   (optional) 1 x M 'cellstr', the spatial variables on which a
%           'polynomial' C may depend;
% - dom:    (optional) M x 2 'double' array, the domains of VARS;
%
% OUTPUTS
% - F:      'polyopvar' object representing the same term as C,
%           expressed as a degree-0 distributed polynomial: a function of
%           the variables of C, or, for a 'dpvar', with the 'dpvar' as
%           coefficient;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - constant2polyopvar
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;

% Allow the input to already be of type 'polyopvar'
if isa(C,'polyopvar')
    F = C;
    return
end

% Determine which independent variables the term depends on
if isa(C,'double')
    pvarname = {};
    C = polynomial(C);
elseif isa(C,'polynomial') || isa(C,'dpvar')
    pvarname = C.varname;
else
    error("Input must be of type 'double', 'polynomial', or 'dpvar'.")
end
if size(C,2)~=1
    error("constant2polyopvar currently only supports column-vector-valued constants.")
end

if isa(C,'dpvar')
    % A decision variable is the coefficient of the degree-0 monomial
    if ~isempty(pvarname)
        error("Decision variables depending on spatial variables are not supported.")
    end
    Cop = C;
else
    % A multiplier operator from R to functions of the variables of C,
    % with their domain taken from vars, dom
    if isempty(pvarname)
        vars = cell(1,0);       dom = zeros(0,2);
    elseif nargin<3
        error("Constants depending on spatial variables require the domain of these variables.")
    else
        is_var = ismember(vars,pvarname);
        if ~all(ismember(pvarname,vars))
            error("The function depends on spatial variables on which the distributed polynomial does not depend.")
        end
        vars = reshape(vars(is_var),1,[]);       dom = dom(is_var,:);
    end
    Cop = mat2tensopvar(C,struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',dom));
end

% Declare a polyopvar representation of the term, as a single degree-0
% monomial (the constant "1") with coefficient Cop
F = polyopvar();
F.varname = cell(0,1);
F.Z = sparse(1,0);
F.C = tensopmat(Cop);
F = initialize(F);

end
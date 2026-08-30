function obj = sopvar2opvar(objSopvar)
% OBJ = SOPVAR2OPVAR(OBJSOPVAR) takes a sopvar object representing a 4-PI
% operator component and returns an opvar object representing the same
% operator.
%
% INPUTS
% - objSopvar:  'sopvar' object representing a 1D PI operator. It can map
%               between different function space, but the input and output
%               domains cannot be both finite- and infinite-dimensional
%               (i.e. L2^n and R^m are supported, but R^m x L2^n is not);
%
% OUTPUTS
% - obj:        'opvar' object representing the same operator as the input;
%
% NOTES
% At this time, the output operator is always expressed in terms of
% variables s1 (primary) and t1 (dummy), independent of the variables
% specified in objSopvar.vars
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - sopvar2opvar
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
% SS, 03/06/2026: Initial coding
% DJ, 05/27/2026: Add support for maps to/from R^n


% Extract the dimension, variables and domain of the operator
dims = objSopvar.dims;
vars = objSopvar.vars;
dom = objSopvar.dom;

% Make sure the operator is 1D
if numel(unique([vars.in,vars.out]))>1
    if numel(unique([vars.in,vars.out]))<=2
        obj = sopvar2opvar2d(objSopvar);
        return
    else
        error("Operator maps between functions of more than one variable.")
    end
end

% Initialize an empty operator
opvar obj;
obj.var1 = polynomial({'s1'});  % use default variables
obj.var2 = polynomial({'t1'});  % should we use s1_dum?

% Distinguish cases of different input/output domains
if isempty(vars.in) && isempty(vars.out)
    % Operator maps R to R
    obj.P = objSopvar.params{1};
elseif isempty(vars.out)
    % Operator maps L2 to R
    obj.I = dom.in;
    Q1 = quadPoly(objSopvar.params{1}, {}, objSopvar.ZR,  dims, {}, {'s1'}, 0);
    Q1 = combine(Q1);
    obj.Q1 = quadPoly.quadPoly2polynomial(Q1);
elseif isempty(vars.in)
    % Operator maps R to L2
    obj.I = dom.out;
    Q2 = quadPoly(objSopvar.params{1}, objSopvar.ZL, {}, dims, {'s1'}, {}, 0);
    Q2 = combine(Q2);
    obj.Q2 = quadPoly.quadPoly2polynomial(Q2);
else
    % Operator maps L2 to L2
    obj.I = objSopvar.dom.in;
    
    R0 = quadPoly(objSopvar.params{1}, objSopvar.ZL, objSopvar.ZR,  dims, {'s1'}, {'s1'}, 0);
    R1 = quadPoly(objSopvar.params{2}, objSopvar.ZL, objSopvar.ZR,  dims, {'s1'}, {'t1'}, 0);   % should we use s1_dum?
    R2 = quadPoly(objSopvar.params{3}, objSopvar.ZL, objSopvar.ZR,  dims, {'s1'}, {'t1'}, 0);
    
    R0 = combine(R0);
    R1 = combine(R1);
    R2 = combine(R2);
    
    obj.R.R0 = combine(quadPoly.quadPoly2polynomial(R0));
    obj.R.R1 = combine(quadPoly.quadPoly2polynomial(R1));
    obj.R.R2 = combine(quadPoly.quadPoly2polynomial(R2)); 
end

end
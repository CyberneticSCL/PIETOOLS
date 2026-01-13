function [P_out] = opvar2ndopvar(P_in,d)
% [P_OUT] = OPVAR2COEFSS(P_IN) takes an opvar object P_IN and returns a
% ndopvar object P_OUT representing the same 3-PI operator.
%
% INPUTS
% - P_in:   mxn 'opvar' object representing a 3-PI operator
% - d:      (optional) scalar 'double' specifiying a desired degree of the
%           monomial basis. Defaults to the maximal degree of the monomials
%           appearing in Pop.R.R0, Pop.R.R1 and Pop.R.R2.
% 
% OUTPUTS
% - P_out:  mxn 'ndopvar' object representing the same 3-PI operator as 
%           P_in, but in the format specified in the 'ndopvar' class file;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opvar2ndopvar
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
% DJ, 01/06/2206: Initial coding

if isa(P_in,'opvar2d') || isa(P_in,'dopvar2d')
    if nargin==1
        P_out = dopvar2d2ndopvar(P_in);
    elseif nargin==2
        P_out = dopvar2d2ndopvar(P_in,d);
    end
    return
elseif ~isa(P_in,'opvar') && ~isa(P_in,'dopvar')
    error("Input  must be of type 'opvar'.")
end
if any(P_in.dim(1,:))
    error("Only operators from L2 to L2 are supported.")
end
if nargin==1
    d = inf;
end

% Extract the relevant information
dim = P_in.dim(2,:);
var1 = P_in.var1;
var2 = P_in.var2;
R0 = dpvar(P_in.R.R0);
R1 = dpvar(P_in.R.R1);
R2 = dpvar(P_in.R.R2);

% Determine which decision variables appear in the operator
dvarname = unique([R0.dvarname; R1.dvarname; R2.dvarname]);

% Determine the maximal monomial degree in both var1 and var2
dmax = max([max(R1.degmat,[],"all"),max(R2.degmat,[],"all")]);
if nargin==1
    d = dmax;
elseif d<dmax
    error("Specified degree is smaller than maximal monomial degree of parameters.")
end

% Get coefficients representing parameters in the quadratic form,
%   R0(s,t) = (Im o Zd(s))^T C0;
%   Ri(s,t) = (Im o Zd(s))^T Ci (In o Zd(t);
C0 = get_quadratic_form(R0,var1,[],dvarname,d);
C1 = get_quadratic_form(R1,var1,var2,dvarname,d);
C2 = get_quadratic_form(R2,var1,var2,dvarname,d);
C_cell = {C0;C1;C2};

% Collect parameters representing the operator in a struct
P_out = ndopvar();
P_out.C = C_cell;
%P_out.dim = dim;
P_out.dom = P_in.I;
P_out.deg = d;
P_out.vars = [var1,var2];
P_out.dvarname = dvarname;

end
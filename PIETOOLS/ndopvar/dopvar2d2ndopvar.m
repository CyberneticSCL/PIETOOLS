function [P_out] = dopvar2d2ndopvar(P_in,d)
% [P_OUT] = DOPVAR2D2NDOPVAR(P_IN,D) takes an opvar2d or dopvar2d object 
% P_IN and returns a ndopvar object P_OUT representing the same 2D PI 
% operator (variable).
%
% INPUTS
% - P_in:   mxn 'opvar2d' or 'dopvard' object representing a 9-PI operator;
% - d:      (optional) scalar 'double' specifiying a desired degree of the
%           monomial basis. Defaults to the maximal degree of the monomials
%           appearing in Pop.R.R1 and Pop.R.R2.
% 
% OUTPUTS
% - P_out:  m x n 'ndopvar' object representing the same 9-PI operator as 
%           P_in, in the format specified in the 'ndopvar' class file;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opvar2d2ndopvar
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

if ~isa(P_in,'opvar2d') && ~isa(P_in,'dopvar2d')
    error("Input  must be of type 'opvar2d'.")
end
if any(any(P_in.dim(1:3,:)))
    error("Only operators from L2 to L2 are supported.")
end
if nargin==1
    d = inf;
end

% Extract the relevant information
dim = P_in.dim(4,:);
var1 = P_in.var1;
var2 = P_in.var2;

% Determine the decision variables and
% maximal monomial degree in var1 and var2
dmax = 0;
P_in = opvar2dopvar2d(P_in);
dvarname = {};
for ii=1:numel(P_in.R22)
    dvarname = [dvarname; P_in.R22{ii}.dvarname];
    dmax = max(dmax,max(max(P_in.R22{ii}.degmat)));
end
if nargin==1
    d = dmax;
elseif d<dmax
    error("Specified degree is smaller than maximal monomial degree of parameters.")
end

% Get coefficients representing the parameters Rij in the qaudratic form,
%   R00(s) = (Im o Zd(s))^T C00;
%   Ri0(s,t1) = (Im o Zd(s))^T Ci0 (In o Zd(t1));
%   R0j(s,t2) = (Im o Zd(s))^T C0j (In o Zd(t2));
%   Rij(s,t) = (Im o Zd(s))^T Cij (In o Zd(t))
% for i,j in {1,2};
C_cell = cell(size(P_in.R22));
C_cell{1} = get_quadratic_form(P_in.R22{1,1},var1,[],dvarname,d);
C_cell{2,1} = get_quadratic_form(P_in.R22{2,1},var1,var2(1),dvarname,d);
C_cell{3,1} = get_quadratic_form(P_in.R22{3,1},var1,var2(1),dvarname,d);
C_cell{1,2} = get_quadratic_form(P_in.R22{1,2},var1,var2(2),dvarname,d);
C_cell{1,3} = get_quadratic_form(P_in.R22{1,3},var1,var2(2),dvarname,d);
for ii=[5,6,8,9]
    C_cell{ii} = get_quadratic_form(P_in.R22{ii},var1,var2,dvarname,d);
end

% Collect parameters representing the operator in a struct
P_out = ndopvar();
P_out.C = C_cell;
%P_out.dim = dim;
P_out.dom = P_in.I;
P_out.deg = [d;d];
P_out.vars = [P_in.var1,P_in.var2];
P_out.dvarname = dvarname;

end
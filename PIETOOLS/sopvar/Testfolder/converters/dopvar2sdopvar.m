function Psop = dopvar2sdopvar(Pop)
% PSOP = DOPVAR2SDOPVAR(POP) takes a dopvar object POP and returns an
% sdopvar object PSOP representing the same decision variable operator
%
% INPUTS
% - Pop:    'dopvar' object representing a single block of a 4-PI operator.
%           That is, only one of Pop.P, Pop.Q1, Pop.Q2, and Pop.R may be
%           nonempty.
%
% OUTPUTS
% - Psop:   'sdopvar' object representing the same operator as Pop;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - dopvar2sdopvar
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
% DJ, 10/05/2026: Initial coding;

% Check that the input is of appropraite class
if isa(Pop,'dopvar2d')
    error("Conversion of 'dopvar2d' objects is currently not supported.")
%    Psop = dopvar2d2sdopvar(Pop);
%    return
elseif isa(Pop,'opvar')
    Pop = opvar2dopvar(Pop);
elseif ~isa(Pop,'dopvar')
    error("Input must be of type 'dopvar'.") 
end
% Make sure the operator maps only to/from one type of space
Pdim = Pop.dim;
if nnz(Pdim(:,1))~=1 || nnz(Pdim(:,2))~=1
    error("Operators mapping between coupled finite- and infinite-dimensional spaces are not supported.")
end

% Determine the variables and domain of the operator
Pdom = Pop.I;
var1 = Pop.var1.varname;    var2 = Pop.var2.varname;
% Initialize empty variables and domain of the output operator
vars = struct;
vars.in = {};           vars.out = {};
dom = struct();
dom.in = zeros(0,2);    dom.out = zeros(0,2);
ZL = cell(1,0);         ZR = cell(1,0);

% Set the parameters, distinguish the four cases depending on
% the input and output spaces of the operator
if Pdim(1,1) && Pdim(1,2)
    % % Pop maps R to R
    dims = Pdim(1,:);
    params = dpvar2sdvar(Pop.P,struct('out',{cell(1,0)},'in',{cell(1,0)}));
    Zd = params.dvarname;
    params = struct('A',{{params.A}},'B',{{params.B}});

elseif Pdim(1,1)
    % % Pop maps L2 to R
    % Extract the kernel
    vars.in = var1;
    [params,~,ZR] = dpvar2sdvar(Pop.Q1,vars);
    Zd = params.dvarname;
    params = struct('A',{{params.A}},'B',{{params.B}});
    % Set the input variables/monomials
    dims = [Pdim(1,1),Pdim(2,2)];
    vars.in = var1;
    dom.in = Pdom;

elseif Pdim(1,2)
    % % Pop maps R to L2
    % Extract the multiplier function
    vars.out = var1;
    [params,ZL,~] = dpvar2sdvar(Pop.Q2,vars);
    Zd = params.dvarname;
    params = struct('A',{{params.A}},'B',{{params.B}});
    % Set the output variables/monomials
    dims = [Pdim(2,1),Pdim(1,2)];
    vars.out = var1;
    dom.out = Pdom;

else
    % % Pop maps L2 to L2
    % Extract the parameters
    vars.in = var2;
    vars.out = var1;
    [A,B,Z,Zd] = dpvars2sdvars({Pop.R.R0;Pop.R.R1;Pop.R.R2},vars);
    params = struct('A',{A},'B',{B});
    % Set the input and output variables/monomials
    dims = Pdim(2,:);
    ZL = Z.out;         ZR = Z.in;
    vars.out = var1;    vars.in = var1;
    dom.out = Pdom;     dom.in = Pdom;
end

% Declare the operator
Psop = sdopvar(params,vars,Zd,ZL,ZR,dom,dims);

end
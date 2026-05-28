function Psop = opvar2sopvar(Pop)
% PSOP = OPVAR2SOPVAR(POP) takes an opvar object POP and returns a sopvar
% object PSOP representing the same operator
%
% INPUTS
% - Pop:    'opvar' object representing a single block of a 4-PI operator.
%           That is, only one of Pop.P, Pop.Q1, Pop.Q2, and Pop.R may be
%           nonempty.
%
% OUTPUTS
% - Psop:   'sopvar' object representing the same operator as Pop;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opvar2sopvar
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
% DJ, 05/27/2026: Initial coding

% Check that the input is of appropraite class
if isa(Pop,'opvar2d')
    Psop = opvar2d2sopvar(Pop);
    return
elseif ~isa(Pop,'opvar')
    error("Input must be of type 'opvar'.")
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
    params = {Pop.P};

elseif Pdim(1,1)
    % % Pop maps L2 to R
    % Extract the kernel
    Qfun = quadPoly.polynomial2quadPoly(Pop.Q1, {}, var1);
    params = Qfun.C;
    params = {full(params)};
    % Set the input variables/monomials
    dims = [Pdim(1,1),Pdim(2,2)];
    ZR = Qfun.Zt;    
    vars.in = var1;
    dom.in = Pdom;

elseif Pdim(1,2)
    % % Pop maps R to L2
    % Extract the multiplier function
    Qfun = quadPoly.polynomial2quadPoly(Pop.Q2, var1, {});
    params = Qfun.C;
    params = {full(params)};
    % Set the output variables/monomials
    dims = [Pdim(2,1),Pdim(1,2)];
    ZL = Qfun.Zs;    
    vars.out = var1;
    dom.out = Pdom;

else
    % % Pop maps L2 to L2
    % Extract the parameters
    R0fun = quadPoly.polynomial2quadPoly(Pop.R.R0, var1, var2);
    R1fun = quadPoly.polynomial2quadPoly(Pop.R.R1, var1, var2);
    R2fun = quadPoly.polynomial2quadPoly(Pop.R.R2, var1, var2);
    [params,Zs,Zt] = get_params({R0fun;R1fun;R2fun});
    % Set the input and output variables/monomials
    dims = Pdim(2,:);
    ZL = Zs;            ZR = Zt;
    vars.out = var1;    vars.in = var1;
    dom.out = Pdom;     dom.in = Pdom;
end

% Declare the operator
Psop = sopvar(params,vars,ZR,ZL,dom,dims);

end



%%
function [C,Zs,Zt] = get_params(Q)
% [C,ZS,ZT] = GET_PARAMS(Q) takes a cell of n quadpoly objects, and
% returns a cell of n coefficient matrices C_OUT and cells ZS and ZT
% of monomials representing the polynomials in Q in terms of the same
% bases.
%
% INPUTS
% - Q:      m x n cell of 'quadpoly' objects, representing matrix-valued
%           polynomials of identical sizes. Must be defined in terms of a
%           single primary variable, and a single dummy variable;
%
% OUTPUTS
% - C:      m x n cell of coefficient matrices, representing the input
%           polynomials in terms of the monomials Zs and Zt
% - Zs:     1 x 1 cell of vectors, specifying a unique basis of monomials
%           in the primary variable in terms of which the input polynomials
%           are expressed using the output coefficients.
% - Zt:     1 x 1 cell of vectors, specifying a unique basis of monomials
%           in the dummy variable in terms of which the input polynomials
%           are expressed using the output coefficients.


% Distinguish case of single polynomial
if isscalar(Q)
    % We can extract the coefficients and monomials directly from Q
    Q = Q{1};
    C = Q.C;        C = {full(C)};
    Zs = Q.Zs;      Zt = Q.Zt;
    return
end

% Extract the coefficients and monomials used in each polynomial
Ccell = cell(size(Q));
Zscell = cell(size(Q));     Ztcell = cell(size(Q));
for i=1:numel(Q)
    Ccell{i} = Q{i}.C;
    Zscell{i} = Q{i}.Zs;   Ztcell{i} = Q{i}.Zt;
end

% Collect the monomials into a big vector
% NOTE: we assume the polynomials are of the same size, and defined in
% terms of a single primary and dummy variable
Zsfull = zeros(0,1);        Ztfull = zeros(0,1);
nZs = zeros(numel(Q),1);    nZt = zeros(numel(Q),1);    % size of each vector
for i=1:numel(Q)
    nZs(i) = size(Zscell{i}{1},1);
    nZt(i) = size(Ztcell{i}{1},1);
    Zsfull = [Zsfull; Zscell{i}{1}];
    Ztfull = [Ztfull; Ztcell{i}{1}];
end
% Declare unique bases of monomials
[Ps,Zs] = uniquerows_integerTable(Zsfull);      % Ps*Zs = Zsfull
[Pt,Zt] = uniquerows_integerTable(Ztfull);
% Augment coefficients with zeros to match new monomial basis
C = Ccell;
[m,n] = size(Q{1});
nnZs = [0;cumsum(nZs)];
nnZt = [0;cumsum(nZt)];
for i=1:numel(C)
    Ps_i = spIkron(m,Ps(nnZs(i)+1:nnZs(i+1),:));
    Pt_i = spIkron(n,Pt(nnZt(i)+1:nnZt(i+1),:));
    C{i} = full(Ps_i'*Ccell{i}*Pt_i);
end
Zs = {Zs};
Zt = {Zt};

end
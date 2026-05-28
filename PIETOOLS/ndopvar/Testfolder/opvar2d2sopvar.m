function Psop = opvar2d2sopvar(Pop)
% PSOP = OPVAR2D2SOPVAR(POP) takes an opvar2d object POP and returns a
% sopvar object PSOP representing the same operator.
%
% INPUTS
% - Pop:    'opvar2d' object representing a single block of a 2D-PI
%           operator. That is, only one of Pop.R00, Pop.R0x, ..., Pop.R22
%           may be nonempty.
%
% OUTPUTS
% - Psop:   'sopvar' object representing the same operator as Pop;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opvar2d2sopvar
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
if isa(Pop,'opvar')
    Psop = opvar2sopvar(Pop);
    return
elseif ~isa(Pop,'opvar2d')
    error("Input must be of type 'opvar2d'.")
end
% Make sure the operator maps only to/from one type of space
Pdim = Pop.dim;
if nnz(Pdim(:,1))~=1 || nnz(Pdim(:,2))~=1
    error("Operators mapping between coupled finite- and infinite-dimensional spaces are not supported.")
end

% Determine the variables and domain of the operator
Pdom = Pop.I;
var1 = pvar2varname(Pop.var1);
var2 = pvar2varname(Pop.var2);
% Initialize empty variables and domain of the output operator
vars = struct;
vars.in = {};           vars.out = {};
dom = struct();
dom.in = zeros(0,2);    dom.out = zeros(0,2);

% Extract non-empty operator from 2D PI structure
ridx = find(Pdim(:,1));
cidx = find(Pdim(:,2));
dims = [Pdim(ridx,1),Pdim(cidx,2)];
Rparam_names = {'R00','R0x','R0y','R02';
                'Rx0','Rxx','Rxy','Rx2';
                'Ry0','Ryx','Ryy','Ry2';
                'R20','R2x','R2y','R22'};
Rname = Rparam_names{ridx,cidx};
R = Pop.(Rname);
if ~isa(R,'cell')
    R = {R};
end
% Make sure dummy variable is used for integration
if ismember(Rname,{'R0x';'R02';'Ry2';'Ryx'})
    for i=1:numel(R)
        R{i} = subs(R{i},Pop.var1(1),Pop.var2(1));
    end
end
if ismember(Rname,{'R0y';'R02';'Rx2';'Rxy'})
    for i=1:numel(R)
        R{i} = subs(R{i},Pop.var1(2),Pop.var2(2));
    end
end

% Determine what variables the operator maps
maps2x = logical(1-mod(ridx,2));
mapsx = logical(1-mod(cidx,2));
maps2y = logical(floor((ridx-1)/2));
mapsy = logical(floor((cidx-1)/2));
vars.in = var1([mapsx;mapsy])';
vars.out = var1([maps2x;maps2y])';
dom.in = Pdom([mapsx;mapsy],:);
dom.out = Pdom([maps2x;maps2y],:);

% Decompose parameters into coefficient matrices and monomials
[params,ZL,ZR] = get_params(R,var1([maps2x;maps2y]),var2([mapsx;mapsy]));
if all(size(params)==[1,3])
    params = params';
end

% Declare the operator
Psop = sopvar(params,vars,ZR,ZL,dom,dims);

end



%%
function [C,Zs,Zt] = get_params(Q,var1,var2)
% [C,ZS,ZT] = GET_PARAMS(Q,VAR1,VAR2) takes a cell of polynomials objects,
% and returns a cell of coefficient matrices representing the polynomials
% in terms of the same bases.
%
% INPUTS
% - Q:      m x n cell of 'polynomial' objects, representing matrix-valued
%           polynomials of identical sizes;
% - var1:   p x 1 cellstr specifying the names of the primary variables in
%           the polynomials;
% - var2:   q x 1 cellstr specifying the names of the dummy variables in
%           the polynomials;
%
% OUTPUTS
% - C:      m x n cell of coefficient matrices, representing the input
%           polynomials in terms of the monomials Zs and Zt
% - Zs:     1 x p cell of vectors, specifying a unique basis of monomials
%           in the primary variable in terms of which the input polynomials
%           are expressed using the output coefficients.
% - Zt:     1 x q cell of vectors, specifying a unique basis of monomials
%           in the dummy variable in terms of which the input polynomials
%           are expressed using the output coefficients.

% Distinguish case of single polynomial
if isscalar(Q)
    % We can extract the coefficients and monomials directly from Q
    Q = quadPoly.polynomial2quadPoly(Q{1}, var1, var2);
    C = Q.C;        C = {C};
    Zs = Q.Zs;      Zt = Q.Zt;
    return
end

% Extract the coefficients and monomials used in each polynomial
Ccell = cell(size(Q));
Zscell = cell(size(Q));     Ztcell = cell(size(Q));
for i=1:numel(Q)
    Q{i} = quadPoly.polynomial2quadPoly(Q{i}, var1, var2);
    C_i = Q{i}.C;
    Ccell{i} = full(C_i);
    Zscell{i} = Q{i}.Zs;   Ztcell{i} = Q{i}.Zt;
end

% Collect the monomials into a big vector
% NOTE: we assume the polynomials are of the same size, and defined in
% terms of a single primary and dummy variable
nvar1 = numel(var1);        nvar2 = numel(var2);
Zsfull = cell(1,nvar1);
Ztfull = cell(1,nvar2);
nZs = zeros(numel(Q),nvar1);    
nZt = zeros(numel(Q),nvar2);    % size of each vector
for i=1:numel(Q)
    for j=1:nvar1
        nZs(i,j) = size(Zscell{i}{j},1);
        Zsfull{j} = [Zsfull{j}; Zscell{i}{j}];
    end
    for j=1:nvar2
        nZt(i,j) = size(Ztcell{i}{j},1);
        Ztfull{j} = [Ztfull{j}; Ztcell{i}{j}];
    end
end
% Declare unique bases of monomials
Zs = cell(1,nvar1);     Ps = 1;
Zt = cell(1,nvar2);     Pt = 1;
for j=1:nvar1
    [Ps_j,Zs_j] = uniquerows_integerTable(Zsfull{j});      % Ps*Zs = Zsfull
    Zs{j} = Zs_j;       
    Ps = kron(Ps,Ps_j);
end
for j=1:nvar2
    [Pt_j,Zt_j] = uniquerows_integerTable(Ztfull{j});
    Zt{j} = Zt_j;       
    Pt = kron(Pt,Pt_j);
end
% Establish which rows of the matrices Ps and Pt correspond to which
% parameters
nnZs = [zeros(1,nvar1);cumsum(nZs,1)];
nnZt = [zeros(1,nvar2);cumsum(nZt,1)];
stride1 = [cumprod(nnZs(end,:),2,'reverse'),1];
stride2 = [cumprod(nnZt(end,:),2,'reverse'),1];
idcs1 = num2cell(ones(size(nZs,1),1));
for i=1:numel(idcs1)
    for j=1:nvar1
        idcs_j = stride1(j+1)*(nnZs(i,j):nnZs(i+1,j)-1)';
        idcs1{i} = idcs1{i}(:)' +idcs_j;
    end
    idcs1{i} = idcs1{i}(:)';
end
idcs2 = num2cell(ones(size(nZt,1),1));
for i=1:numel(idcs2)
    for j=1:nvar2
        idcs_j = stride2(j+1)*(nnZt(i,j):nnZt(i+1,j)-1)';
        idcs2{i} = idcs2{i}(:)' +idcs_j;
    end
    idcs2{i} = idcs2{i}(:)';
end


% Augment coefficients with zeros to match new monomial basis
C = Ccell;
[m,n] = size(Q{1}); % assume same size for all polynomials
for i=1:numel(C)
    Ps_i = spIkron(m,Ps(idcs1{i},:));
    Pt_i = spIkron(n,Pt(idcs2{i},:));
    C{i} = full(Ps_i'*Ccell{i}*Pt_i);
end

end
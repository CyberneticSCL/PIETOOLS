function P = PIETOOLS_random_poly_generator(nr,nc,varname,nZ,dmax)
% P = PIETOOLS_random_poly_generator(nr,nc,varname,dmax) generates a
% ``random'' nr x nc 'polynomial' class objects in variables "varname",
% built on a basis of at most nZ monomials of degree at most dmax.
%
% INPUT
% - nr:      1x1 integer specifying the row-dimension of the desired object.
% - nc:      1x1 integer specifying the column-dimension of the desired object.
% - varname: cellstr type object specifying the names of the variables on
%            which the desired polynomial can depend.
% - nZ:      1x1 integer specifying the maximal number of distinct 
%            monomials that may appear in the polynomial.
% - dmax:    1x1 integer specifying the maximally allowed combined degree
%            of each monomial that appears in the polynomial.
%
% OUTPUT:
% - P:       nr x nc polynomial class object with P.varname = varname. The
%            coefficients of this monomial are randomly generated,
%            integer-valued. The monomials that appear are randomly
%            selected from all possible monomials in the proposed variables
%            up to maximal degree dmax, including no more than nZ distinct
%            monomials.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_random_poly_generator
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 08/23/2022
%

% % % Process the inputs
if nargin<=2
    error('Insufficient input arguments.')
elseif nargin==2
    if isa(nc,'polynomial') || iscellstr(nc) || ischar(nc)
        varname = nc;
        if numel(nr)==1
            nc = nr;
        else
            nc = nr(2);
            nr = nr(1);
        end
        nZ = 3;
        dmax = 3;
    else
        error('No variable names have been provided.')
    end
elseif nargin==3
    if isa(nc,'polynomial') || iscellstr(nc) || ischar(nc)
        nZ = varname;
        varname = nc;
        if numel(nr)==1
            nc = nr;
        else
            nc = nr(2);
            nr = nr(1);
        end
    end
    dmax = nZ;
elseif nargin==4
    if isa(nc,'polynomial') || iscellstr(nc) || ischar(nc)
        nZ = varname;
        varname = nc;
        if numel(nr)==1
            nc = nr;
        else
            nc = nr(2);
            nr = nr(1);
        end
    end
    dmax = nZ;
end

% % % Error check the inputs
% Check the dimensions
if numel(nr)~=1 || numel(nc)~=1 || any([nr,nc]<0)
    error('Matrix dimensions nr and nc should be specified as separate, 1x1, nonnegative integers.')
end

% Convert varname to cellstr.
if isempty(varname)
    varname = cell(0,1);
elseif ischar(varname)
    varname = {varname};
elseif ispvar(varname)
    varname = varname.varname;
else
    error('Variable names should be specified as pvar or cellstr type object.');
end

nvars = length(varname);
matdim = [nr nc];

% Establish degrees for the variables in each possible monomial up to
% maximal degree dmax;
degmat = monomials(nvars,dmax);
if nZ>size(degmat,1)
    warning(['The desired size of the monomial basis exceeds the maximal amount of monomials possible in the specified variables up to the specified degree.'])
    nZ = size(degmat,1);
else
    indcs = randperm(size(degmat,1));
    degmat = degmat(indcs(1:nZ),:);
end

% Generate random integer coefficients up to values dmax+1;
nm = nr*nc;
coeff = randi([0,dmax+1],nZ,nm);

% Build the polynomial
P = polynomial(coeff,degmat,varname,matdim);

end
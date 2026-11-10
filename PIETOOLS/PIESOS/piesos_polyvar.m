function [prog,p] = piesos_polyvar(prog,Zx,opts)
% [PROG,P] = PIESOS_POLYVAR(PROG,ZX,OPTS) adds an indefinite distributed
% polynomial functional variable P in distributed monomials ZX to the
% PIESOS program structure PROG.
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program (see `piesos_program`);
% - Zx:     nZ x 1 'polyopvar' object representing a basis of nZ
%           distributed monomials in terms of which the polynomial
%           variable "p" is to be defined. A functional acts on each
%           monomial i, defined on the function space of the output of
%           row i of Zx, over the spatial variables on which this row
%           depends;
% - opts:   (optional) struct specifying options for the operator C
%           defining the polynomial, p = C*Zx, with fields
%               'intmat':  OPTION 1: nZ x 1 cell, specifying for each of
%                            the monomials the intmat of the functional
%                            acting on this operator in one of the two
%                            formats specified in `piesos_intopvar`;
%                          OPTION 2: 'char' object, one of 'full',
%                            'no_mult', 'lower' or 'upper', specifying
%                            what type of intmat to use for all of the
%                            distributed monomials. See the
%                            `piesos_intopvar` function;
%                          Here, 'no_mult' is passed on as 'nomult'.
%                          Defaults to 'full'. A cell intmat for a
%                          monomial must have an element for each of the
%                          N variables in Zx.vars.out, where the
%                          elements for variables on which the monomial
%                          does not depend are ignored;
%               'degs':    degrees of the monomials in the kernels of the
%                          functionals, as `degs` in `piesos_intopvar`:
%                          a struct with field 'max', the maximal
%                          degrees, or 'full', the full arrays of
%                          monomial degrees. A cell or array is read as
%                          struct('max',degs). It can be specified as
%                          - nZ x 1 cell, with element i the degrees for
%                            monomial i;
%                          - a single struct, scalar or 1 x N array, for
%                            all monomials, where a scalar is the
%                            maximal cumulative degree in the dummy
%                            variables of each of the N variables in
%                            Zx.vars.out, and a 1 x N array the maximal
%                            cumulative degree for each variable;
%                          - nZ x 1 array, with a maximal cumulative
%                            degree for each monomial.
%                          As for the intmat, the degrees for a monomial
%                          are specified for each of the N variables: a
%                          1 x N cell for 'max', where the elements for
%                          variables on which the monomial does not
%                          depend are ignored; for 'full' a 1 x N cell
%                          of the monomials of each variable, or an
%                          ntrms(1) x ... x ntrms(N) cell of the terms,
%                          with singleton dimensions for these
%                          variables, where such a variable has an array
%                          without columns (e.g. []) in a 1 x N cell.
%                          Defaults to 0 (constant kernels);
%
% OUTPUTS
% - prog:   'struct' specifying the same PIESOS program as the input, but
%           with the polynomial variable p added;
% - p:      scalar 'polyopvar' object representing a polynomial functional
%           parameterized by decision variables. It is defined as
%               p = Cop*Zx,    Cop = [Kop_1,...,Kop_nZ]
%           where Kop_i is the functional acting on monomial i, declared by
%           `piesos_intopvar`;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_polyvar
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
% DJ, 10/15/2026: Initial coding;


% Check the inputs
if ~isa(Zx,'polyopvar')
    error("Basis of distributed monomials must be specified as 'polyopvar' object.")
end
if ~Zx.chkval
    Zx = initialize(Zx);        % fill in the function spaces of the monomials
end

% Check the options
if nargin<=2
    opts = struct();
end
if ~isa(opts,'struct')
    error("Options must be specified as a 'struct'.")
end
nZ = size(Zx.C.ops,1);
N = numel(Zx.vars.out);

% Take the intmat of each of the nZ monomials
if ~isfield(opts,'intmat')
    intmat = repmat({'full'},nZ,1);
elseif ischar(opts.intmat) || (isstring(opts.intmat) && isscalar(opts.intmat))
    intmat = repmat({char(opts.intmat)},nZ,1);
elseif isa(opts.intmat,'cell') && numel(opts.intmat)==nZ
    intmat = opts.intmat(:);
else
    error("Option 'intmat' must be specified as 'char' or as nZ x 1 cell, for nZ distributed monomials.")
end
% Take the degrees of each of the nZ monomials
if ~isfield(opts,'degs')
    % Assume 0 degrees by default
    degs = repmat({0},nZ,1);
elseif isstruct(opts.degs) || (isnumeric(opts.degs) && (isscalar(opts.degs) || (isrow(opts.degs) && numel(opts.degs)==N)))
    % The same degrees for all monomials
    degs = repmat({opts.degs},nZ,1);
elseif isnumeric(opts.degs) && iscolumn(opts.degs) && numel(opts.degs)==nZ
    % A cumulative degree for each monomial
    degs = num2cell(opts.degs);
elseif isa(opts.degs,'cell') && numel(opts.degs)==nZ
    degs = opts.degs(:);
else
    error("Degrees must be specified as struct, scalar, 1 x N array, or nZ x 1 array or cell, for N spatial variables and nZ distributed monomials.")
end

% Declare an operator variable acting on the distributed monomials
Kcell = cell(1,nZ);
for i=1:size(Zx.C.ops,1)
    % Determine the function space in which monomial i exists
    dim_i = Zx.dims.out(i);
    deps_i = Zx.depmat.out(i,:);
    % Declare an operator variable acting on this function space
    intmat_i = intmat{i};
    if ischar(intmat_i) || isstring(intmat_i)
        if strcmp(intmat_i,'no_mult')
            intmat_i = 'nomult';
        end
    elseif isa(intmat_i,'cell') && numel(intmat_i)==N
        % Variables with deps_i=0 are dropped in piesos_intopvar: no columns
        intmat_i = reshape(intmat_i,1,N);
        intmat_i(deps_i==0) = {zeros(1,0)};
    else
        error("The intmat of monomial %d must be specified as 'char' or 1 x N cell, for N spatial variables.",i)
    end
    degs_i = degs{i};
    if isnumeric(degs_i) && ~isscalar(degs_i)
        % Maximal cumulative degree for each variable, as the cell of 'max'
        degs_i = num2cell(reshape(degs_i,1,[]));
    end
    [prog,Kop_i] = piesos_intopvar(prog,[1,dim_i],Zx.vars.out,Zx.dom.out,deps_i,intmat_i,degs_i);
    Kcell{i} = Kop_i;
end
% Collect the functionals in a row of operators, acting on the rows of Zx
Cop = tensopmat(Kcell);
% Compose with the operators in Zx.C, from the tensor monomials to the rows of
% Zx: the identity for a degree-1 monomial, but for a product of states
% the intop*tensopvar composition, which restricts the kernel to the diagonal
% (a delta term)
p = Cop*Zx;



end

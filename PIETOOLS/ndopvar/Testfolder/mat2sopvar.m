function Pop = mat2sopvar(Pmat,vars,dom)
% POP = MAT2SOPVAR(PMAT,VARS,DOM) takes a matrix PMAT and returns a
% 'sopvar' object POP representing the associated multiplier operator
% between the function spaces defined by DOM.
%
% INPUTS
% - Pmat:   m x n 'double' or 'polynomial' array, representing a multiplier
%           operator
% - vars:   struct with fields
%     in:   1 x N cellstr array specifying the names of the input variables
%           to the operator;
%    out:   1 x M cellstr array specifying the names of the output
%           variables to the operator;
%
% - dom:    struct with fields
%    in:    N x 2 'double' array specifying the intervals on which the
%           the input variables is defined;
%   out:    M x 2 'double' array specifying the intervals on which the
%           output variables are defined;
%
% OUTPUTS
% - Pop:    m x n 'sopvar' object representing the multiplier operator
%           defined by Pmat;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mat2sopvar
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
% DJ, 06/12/2026: Initial coding

% Check that the matrix is properly specified
if isa(Pmat,'double')
    Pmat = polynomial(Pmat);
elseif ~isa(Pmat,'polynomial')
    error("Matrix must be specified as object of type 'double' or 'polynomial'.")
end

% Check that the variables are properly specified
if ~isa(vars,'struct') || ~isfield(vars,'in') || ~isfield(vars,'out')
    error("Variables must be specified as struct with fields 'in' and 'out'.")
end
if ~iscellstr(vars.in) || ~iscellstr(vars.out)
    error("Input and output variable names should be specified as 'cellstr' objects.")
end
vars.in = reshape(vars.in,1,[]);        N = numel(vars.in);
vars.out = reshape(vars.out,1,[]);      M = numel(vars.out);

% Check that the domains are properly specified
if ~isa(dom,'struct') || ~isfield(dom,'in') || ~isfield(dom,'out')
    error("Domains must be specified as struct with fields 'in' and 'out'.")
end
if ~isnumeric(dom.in) || ~all(size(dom.in)==[N,2])
    error("Number of input domains should match number of input variables.")
end
if ~isnumeric(dom.out) || ~all(size(dom.out)==[M,2])
    error("Number of output domains should match number of output variables.")
end

% Extract the coefficients defining the multiplier operator
Pmat = quadPoly.polynomial2quadPoly(Pmat,vars.out,cell(1,0));
Cmat = Pmat.C;
ZL = Pmat.Zs;       ZR = num2cell(zeros(1,N));

% Declare zero parameters for PI operator from vars.in to vars.out
vars_S3 = intersect(vars.in,vars.out);
N_S3 = numel(vars_S3);
params = repmat({spalloc(size(Cmat,1),size(Cmat,2),0)},[3*ones(1,N_S3),1]);

% Set multiplier term
params{1} = Cmat;

% Declare the operator
Pop = sopvar(params,vars,ZR,ZL,dom,size(Pmat));

end
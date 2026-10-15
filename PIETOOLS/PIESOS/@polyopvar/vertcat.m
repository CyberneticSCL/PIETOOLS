function F = vertcat(varargin)
% F = VERTCAT(VARARGIN) takes the vertical concatenation of distributed
% polynomials listed in VARARGIN, e.g.
%   F = [varargin{1}; varargin{2}; ...; varargin{end}];
%
% INPUTS
% - varargin:   'polyopvar' objects representing distributed polynomials,
%               or numeric arrays, which are converted to 'polyopvar'
%               objects using constant2polyopvar. Empty inputs are
%               ignored;
%
% OUTPUTS
% - F:          'polyopvar' object representing the vertical
%               concatenation of the input functions;
%
% NOTES
% - The input function spaces of the coefficients defining the
%   polynomials (C.dims.in, C.vars.in and C.depmat.in) must match.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - vertcat
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

% Return an empty 'polyopvar' object if no inputs are specified
if nargin==0
    F = polyopvar();
    return
end

% Check that the first input is appropriate
F = varargin{1};
if isempty(F)
    F = vertcat(varargin{2:end});
    return
end
if isnumeric(F)
    F = constant2polyopvar(F);
elseif ~isa(F,'polyopvar')
    error("Concatenation is only supported with objects of type 'polyopvar'.")
end
if ~F.chkval
    F = initialize(F);
end
F.chkval = true;

% Concatenate with the remaining functions
for i=2:nargin
    % Extract the ith function
    B = varargin{i};
    if isempty(B)
        continue
    elseif isnumeric(B)
        B = constant2polyopvar(B);
    elseif ~isa(B,'polyopvar')
        error("Concatenation is only supported with objects of type 'polyopvar'.")
    end
    if ~B.chkval
        B = initialize(B);
    end
    % Express the operators in terms of the same basis
    [F,B] = common_basis(F,B);
    % Make sure the output domains match
    if ~isequal(B.C.dims.in,F.C.dims.in)...
        || ~isequal(B.C.vars.in,F.C.vars.in)...
            || ~isequal(B.C.depmat.in,F.C.depmat.in)
        error("Input function spaces of the coefficients defining the polynomials must match for vertical concatenation.")
    end
    % Concatenate the coefficient operators defining each polynomial
    F.C = [F.C; B.C];
    F.chkval = true;
    % Set the output domains of the function
    F.dims.out = [F.dims.out; B.dims.out];
    F.depmat.out = [F.depmat.out; B.depmat.out];
    F.chkval = true;
end

end
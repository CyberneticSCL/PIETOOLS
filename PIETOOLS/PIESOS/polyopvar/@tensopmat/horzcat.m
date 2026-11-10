function Pop = horzcat(varargin)
% POP = HORZCAT(VARARGIN) takes the horizontal concatenation of tensor-PI 
% operators listed in VARARGIN, e.g.
%   POP = [varargin{1}, varargin{2}, ..., varargin{end}];
%
% INPUTS
% - varargin:   'tensopmat' objects representing arrays of tensor-PI
%               operators;
%
% OUTPUTS
% - Pop:        'tensopmat' object representing the horizontal
%               concatenation of the input operators;
%
% NOTES
% - Inputs of type 'tensopvar', 'intop' or 'sopvar' are converted to
%   'tensopmat', and empty inputs are skipped.
% - The output function spaces of the operators must match.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - horzcat
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

% Check that the first input is appropriate
Pop = varargin{1};
if isempty(Pop)
    args = varargin(2:end);
    if isempty(args)
        return
    end
    Pop = horzcat(args{:});
    return
end
if isa(Pop,'sopvar') || isa(Pop,'tensopvar') || isa(Pop,'intop')
    Pop = tensopmat(Pop);
elseif ~isa(Pop,'tensopmat')
    error("Concatenation is only supported with objects of type 'tensopmat', 'tensopvar', or 'intop'.")
end
if ~Pop.chkval && ~isvalid(Pop)
    error("Concatenation is not supported: One of the 'tensopmat' objects is not properly declared.")
end
Pop.chkval = true;

% Concatenate with the remaining operators
for i=2:nargin
    % Extract the ith operator
    Bop = varargin{i};
    if isempty(Bop)
        continue
    end
    if isa(Bop,'sopvar') || isa(Bop,'tensopvar') || isa(Bop,'intop')
        Bop = tensopmat(Bop);
    elseif ~isa(Bop,'tensopmat')
        error("Concatenation is only supported with objects of type 'tensopmat', 'tensopvar', or 'intop'.")
    end
    if ~Bop.chkval && ~isvalid(Bop)
        error("Concatenation is not supported: One of the 'tensopmat' objects is not properly declared.")
    end
    Bop.chkval = true;
    % Express the operators in terms of the same variables
    [Pop,Bop] = common_vars(Pop,Bop);
    % Rows of blocks of type double/polynomial/dpvar map to the variables 
    % appearing in them
    [Pop,Bop] = merge_out_spaces(Pop,Bop);
    % Make sure the output domains match
    if ~isequal(Bop.dims.out,Pop.dims.out)...
        || ~isequal(Bop.vars.out,Pop.vars.out)...
            || ~isequal(Bop.depmat.out,Pop.depmat.out)
        error("Output function spaces of the operators must match for horizontal concatenation.")
    end
    % Concatenate the operators
    Pop.ops = [Pop.ops, Bop.ops];
    % Set the input domains of the operator
    Pop.dims.in = [Pop.dims.in; Bop.dims.in];
    Pop.depmat.in = [Pop.depmat.in; Bop.depmat.in];
    Pop.chkval = true;
end

end
function [PDE_out] = horzcat(varargin)
% PDE_OUT = HORZCAT(PDE_1,PDE_2) horizontally concatenates the PDEs defined
% by PDE_1 and PDE_2. This comes down to adding the free terms in PDE_2 to
% either the equations specified in PDE_1, or the free terms specified in
% PDE_1.
% 
% INPUT
% - PDE_1, PDE_2, PDE_3, ...:   
%                   'pde_struct' objects. The first object can represented
%                   completed PDE equations, to which to add some terms.
%                   All other objects must correspond to free terms to be
%                   used to construct PDEs.
%
% OUTPUT
% - PDE_out:        'pde_struct' object containing the same equations or 
%                   terms as in PDE_1, but now with the terms added from
%                   all other input PDE structures.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 06/23/2024


% % % Process the inputs
% % If more than two inputs are provided, just vertcat repeatedly.
if nargin==1
    % Vertcat of a single object is just the object;
    PDE_out = varargin{1};
    return
elseif nargin==2
    % Vertcat of two objects we process below;
    PDE_1 = varargin{1};    PDE_2 = varargin{2};
else
    % Vertcat of more than two objects we just keep repeating;
    PDE_1 = varargin{1};    PDE_2 = varargin{2};
    PDE_rem = varargin{3:end};
end

% % Check that the inputs are appropriate.
if ~isa(PDE_1,'pde_struct')
    error("Concatenation of 'pde_struct' objects with non-pde struct objects is not supported.")
end

% % Perform the concatenation: just use plus to add the terms.
PDE_out = plus(PDE_1,PDE_2);
if nargin>=3
    PDE_out = horzcat(PDE_out,PDE_rem);
end

end
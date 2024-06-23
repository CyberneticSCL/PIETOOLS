function PDE_out = uplus(PDE_in)
% PDE_OUT = UPLUS(PDE_IN) takes the unary plus (+PDE_in) of the input PDE.
% This just returns the same PDE structure...
%
% INPUT
% - PDE_in:     'pde_struct' object representing representing PDE and
%               output equations, or loose terms (PDE_in.free) yet to be
%               used to construct an actual PDE;
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing the same equations or
%               terms as the input. Yes, really.
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

PDE_out = PDE_in;

end
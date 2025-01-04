function PDE_out = uminus(PDE_in)
% PDE_OUT = UMINUS(PDE_in) takes a PDE structure "PDE_in" and returns a
% structure "PDE_out" representing the same PDE and same free terms but
% with sign-swapped coefficients.
%
% INPUT
% - PDE_in:     'pde_struct' object representing representing PDE and
%               output equations, or loose terms (PDE_in.free) yet to be
%               used to construct an actual PDE;
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing the same equations or
%               terms as the input, but now with all coefficients
%               sign-swapped.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024 PIETOOLS Team
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
% DJ, 06/23/2024: Initial coding;
% DJ, 01/03/2025: Update to assume a loose PDE variable is specified as a
%                   single free term, see also update to "pde_var";

% % Process the input.
% Check that the input is indeed of suitable type.
if ~isa(PDE_in,'pde_struct')
    error("This function should only be called for 'pde_struct'-type objects...")
end

% % Change the sign of all the coefficients of all the terms in all of the
% % equations of each equation type.
PDE_out = PDE_in;
eq_types = {'free';'x';'y';'z'};
for kk=1:numel(eq_types)
    eq_type = eq_types{kk};
    for ii=1:numel(PDE_out.(eq_type))
        if ~isfield(PDE_out.(eq_type){ii},'term') || isempty(PDE_out.(eq_type){ii}.term)
            % There are no terms of which to change the sign
            continue
        else
            for jj=1:numel(PDE_out.(eq_type){ii}.term)
                if isfield(PDE_out.(eq_type){ii}.term{jj},'C')
                    PDE_out.(eq_type){ii}.term{jj}.C = -PDE_in.(eq_type){ii}.term{jj}.C;
                else
                    PDE_out.(eq_type){ii}.term{jj}.C = -1;
                end
            end
        end
    end
end

end
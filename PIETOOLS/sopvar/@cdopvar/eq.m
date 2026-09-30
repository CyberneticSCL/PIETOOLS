function tf = eq(A,B,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TF = EQ(A,B,TOL) or A==B tests whether two containers represent the same
% operator; A==0 tests whether A is zero.
%
% INPUTS
% - A, B:   'cdopvar' or 'copvar' objects, single 'sopvar'/'sdopvar'
%           blocks, or the double 0;
% - tol:    (optional) absolute coefficient tolerance, default 1e-14 as for
%           the blocks;
%
% OUTPUTS
% - tf:     logical scalar; true iff the same ordered sequence of component % MMP, 09/29/2026
%           spaces (by name and domain) on each side and equal kernels      % MMP, 09/29/2026
%           componentwise, within tol, whatever the two block partitions,   % MMP, 09/29/2026
%           as affine functions of the decision variables, matched by       % MMP, 09/29/2026
%           name. Different spaces give false, not an error.                % MMP, 09/29/2026
%
% NOTES
% The logic is in 'eq_copvar', shared with 'copvar'. The two containers may
% hold different decision variable lists: the block 'eq' aligns them.
%
% See also EQ_COPVAR, MINUS, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eq(cdopvar)
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
% Initial coding MMP, 09/25/2026
% MMP, 09/29/2026: OUTPUTS updated for 'eq_copvar', which now compares one
%                  operator stored in two block partitions as equal. No code
%                  change here. The OUTPUTS entry for tf was:
% - tf:     logical scalar; true iff same spaces, dimensions and domains,   % MMP, 09/29/2026 (was)
%           and every block pair equal within tol as affine functions of    % MMP, 09/29/2026 (was)
%           the decision variables, matched by name. Different spaces give  % MMP, 09/29/2026 (was)
%           false, not an error.                                            % MMP, 09/29/2026 (was)

if nargin<3,    tol = [];   end
tf = eq_copvar(A,B,tol);

end

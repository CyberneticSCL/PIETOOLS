function out_pie = ctranspose(in_pie)
% OUT_PIE = CTRANSPOSE(IN_PIE) constructs and returns the dual of a given
% PIE system. Given a pie_struct with parameters of the PIE
%   T \dot x = A x + Bw w,      z = Cz x + Dzw w
% the output is the pie_struct with parameters corresponding to the PIE
%   T' \dot x = A' x + Cz' w,   z = Bw' x + Dzw' w
% i.e. out_pie.T = in_pie.T', out_pie.A = in_pie.A',
% out_pie.Cz = in_pie.Bw', out_pie.Bw = in_pie.Cz', out_pie.Dzw =
% in_pie.Dzw'.
%
% INPUT
% - in_pie:     pie_struct class object with Tw=0 and Tu=0, and without
%               algebraic constraint (K, Kw, and Ku without any rows);
%
% OUTPUT
% - out_pie:    pie_struct class object that represents the dual of the
%               in_pie.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ctranspose
%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding SS  - 9_7_2023
% DJ, 10/29/2026: Update to new pie_struct structure (copvar operators, and
%                   no dual for PIEs with an algebraic constraint);


% The operators are copvar objects, so the PIE must be initialized, and the
% function spaces of the exogenous input w and regulated output z are
% exchanged. 
if ~in_pie.is_initialized
    in_pie = initialize(in_pie);
end
% The operators Tw, Dzu and Dyw have no counterpart in the dual.
if ~(in_pie.Tw==0) || ~(in_pie.Tu==0)
    error('Dual PIE cannot be found for systems with disturbances at the boundary');
end
% The dual is not defined for PIEs with an algebraic constraint
% Ku*u + Kw*w + K*x = 0 on the state.
if size(in_pie.K,1)>0
    error('The dual PIE is not supported for PIEs with an algebraic constraint (K, Kw, Ku).');
end

out_pie = in_pie;
out_pie.A = in_pie.A';
out_pie.T = in_pie.T';
out_pie.Cz = in_pie.Bw';
out_pie.Bw = in_pie.Cz';
out_pie.Dzw = in_pie.Dzw';
out_pie.Tw = [];
out_pie.Dzu = [];
out_pie.Dyw = [];
out_pie.K = [];     out_pie.Kw = [];    out_pie.Ku = [];
out_pie.space.w = in_pie.space.z;
out_pie.space.z = in_pie.space.w;
% The dual PIE no longer corresponds to the original PDE.
default_pie = pie_struct();
out_pie.space_pde = default_pie.space_pde;
out_pie.pie2pde = default_pie.pie2pde;
out_pie.pde2pie = default_pie.pde2pie;
out_pie = initialize(out_pie);
end
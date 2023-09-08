function out_pie = ctranspose(in_pie)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_pie = ctranspose(in_pie) constructs and returns the dual of a given
% PIE system. For e.g.,
% Given a pie_struct with parameters of the PIE
% T \dot x = A x+ B1 w, z = C1 x + D11 w
% the output is the pie_struct with parameters corresponding to the PIE
% T' \dot x = A' x+ C1' w, z = B1' x + D11' w
% i.e., out_pie.T = in_pie.T', out_pie.A = in_pie.A', out_pie.C1 = in_pie.B1'
%        out_pie.B1 = in_pie.C1', out_pie.D11 = in_pie.D11'
% 
% INPUT
% in_pie: pie_struct class object with Tw=0 and Tu=0
% 
% OUTPUT
% out_pie: pie_struct class object that represents the dual of the in_pie 
% 
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2023  M. Peet, S. Shivakumar, D. Jagt
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

if ~(in_pie.Tw==0)||~(in_pie.Tu==0)
    error('Dual PIE cannot be found for systems with disturbances at the boundary');
end

out_pie=in_pie;
out_pie.A= in_pie.A';
out_pie.T= in_pie.T';
out_pie.C1= in_pie.B1';
out_pie.B1= in_pie.C1';
out_pie.D11= in_pie.D11';
end
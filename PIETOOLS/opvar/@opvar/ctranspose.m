function [Pt] = ctranspose(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P] = ctranspose(P) transposes an operator P: R^p x L2^q to R^m x L2^n
% Date: 6/13/19
% Version: 1.0
% 
% INPUT
% P: opvar class object
% 
% OUTPUT
% Padj: transpose of the input opvar with the same matlab structure as
% 
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - ctranspose
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 7_26_2019
%
if ~isa(P,'opvar')
    error('Input must be an opvar variable.');
end

opvar Pt; Pt.I = P.I; Pt.var1= P.var1; Pt.var2 = P.var2;

Pt.P = P.P';
Pt.Q2 = P.Q1';
Pt.Q1 = P.Q2';
Pt.R.R0 = P.R.R0';
Pt.R.R2 = var_swap(P.R.R1',P.var1, P.var2);
Pt.R.R1 = var_swap(P.R.R2',P.var1, P.var2);
end
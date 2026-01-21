function C = plus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = plus(A,B) adds two sopvar operators A,B: L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% B: sopvar class object L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% Note: A,B must have identical dimensions and domains
%
% OUTPUT
% C = A+B:  L_2^q[Sj,...,Sm] to L_2^p[Si,...Sn]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - plus
%
% Copyright (C)2026  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 1_16_2026
%

% initialize added output class
C = B;

% Error handling: Checks to ensure A and B are compatible
if any(A.dims~=B.dims)
    error('Summands A and B have different dimensions');
end
if any(~strcmp(A.vars_in,B.vars_in)) ||any(~strcmp(A.vars_out,B.vars_out))
    error('Summands A and B map between different spaces');
end
if any(any(A.dom_in~=B.dom_in)) || any(any(A.dom_out~=B.dom_out))
    error('input or output variables in summands A and B have different domains');
end

if numel(A.params)~=numel(B.params)
    error('number of terms in summands is not equal -- one of them is probably malformed');
end
% dimensions of params in A and B, vars_in and vars_out

for i=1:numel(A.params)  % linear indexing of multi-dimensional cell array
    C.params{i} = A.params{i}+C.params{i};  % adding quadpoly objects. % Note: chatgpt says this approach is faster than cellfun
end
end

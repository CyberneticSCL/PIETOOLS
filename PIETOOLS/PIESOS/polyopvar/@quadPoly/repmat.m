function P_repmat = repmat(P1, N, M)
% F = repmat(varargin) generates repmat of quadpoly objects using
% varargs as the blocks on the diagonal
% 
% INPUTS
% - varargin:   'quadpoly' class objects.
% 
% OUTPUTS
% - P_repmat:      'quadpoly' object with each of the parameters definined as the
%               repmat of quadpoly 
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  PIETOOLS Team
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
% AT, 02/18/2026: Initial coding;

if ~isa(P1, 'quadPoly')
    error('quadPoly:repmat:badType', 'repmat supports quadPoly only.');
end

if nargin==1
    P_repmat=P1;
    return 
end



if nargin==2
    M = N;
end
 


% extract monomials and coefficients.
Zt = P1.Zt;
Zs = P1.Zs;
nt = P1.nt;
ns = P1.ns;
C = P1.C;

C_repmat = repmat(C, N, M); 

P_repmat = quadPoly(C_repmat, Zs, Zt, P1.dim.*[N, M], ns, nt); 

end
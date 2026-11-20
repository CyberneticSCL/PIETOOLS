function Pop = int2copvar(Kmat,vars,dom,int_dom)
% POP = INT2COPVAR(KMAT,VARS,DOM,INT_DOM) takes a kernel KMAT and returns a
% 'copvar' object POP representing the associated integral operator between
% the function spaces defined by DOM,
%   (Pop*x)(s) = int_{INT_DOM} Kmat(s,t) x(t) dt,
% with s the variables vars.out, and t the variables vars.in.
%
% INPUTS
% - Kmat:   m x n 'double', 'polynomial' or 'quadPoly' array, the kernel
%           K(s,t); in a 'polynomial', the dummy variable t of a variable s
%           that is also in vars.in is named s_dum;
% - vars:   struct with fields 'in' (1 x N cellstr, the variables t) and
%           'out' (1 x M cellstr, the variables s), or one cellstr for both;
% - dom:    struct with fields 'in' (N x 2 'double') and 'out' (M x 2
%           'double'), the intervals of the variables, or one N x 2 array
%           for both;
% - int_dom: (optional) N x 2 'double', 'polynomial' or 'cell' array, the
%           limits [L,U] of integration in each variable of vars.in: the
%           limits of its domain, or the variable itself if it is also in
%           vars.out. Defaults to the full domain;
%
% OUTPUTS
% - Pop:    1 x 1 'copvar' object holding the 'sopvar' integral operator
%           defined by Kmat;
%
% NOTES
% - The operator has a single input space, with the variables vars.in, and
%   a single output space, with the variables vars.out;
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - int2copvar
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
% DJ, 10/30/2026: Initial coding;

% Default limits: the full domain
if nargin<=3
    int_dom = [];
end

% Declare the operator as sopvar, and use it as the block of the container
Pop = copvar({int2sopvar(Kmat,vars,dom,int_dom)});

end

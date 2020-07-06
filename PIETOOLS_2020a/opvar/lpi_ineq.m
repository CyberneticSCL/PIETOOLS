function sos = lpi_ineq(sos,P, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sos = lpi_ineq(sos,P, options) sets up inequality constraints such that P>=0
% 
% INPUT
%   prog: SOS program to modify.
%   P: PI opvar variable
%   -Optional INPUTS
%   options: determines if the operator needs to be pure or to include psatz
%   term. If options.pure=1, then R_0 term is excluded. If options.psatz=1,
%   then the psatz term is included
% 
% OUTPUT 
%   sos: SOS program
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpi_ineq
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
% Initial coding MMP, SS  - 9_13_2019
%

switch nargin
    case 1
        error(['Not enough inputs!'])
    case 2
        options.psatz=0;
        options.pure=0;
        options.sep =0;
    case 3
        if ~isfield(options,'psatz')
            options.psatz=0;
        end
        if ~isfield(options,'pure')
            options.pure=0;
        end
        if ~isfield(options,'sep')
            options.sep=0;
        end
end

dim = P.dim;
if dim(:,1)~=dim(:,2)
    error('Non-symmetric Operators cannot be sign definite. Unable to set the inequality');
end

nx1 = dim(1,1);
nx2 = dim(2,1);
X = P.I;


d2 = degbalance(P);
if options.pure == 1
    options2.exclude= [0 1 0 0];
    options3.exclude = [0 1 0 0]; 
else
    options2.exclude= [0 0 0 0];
    options3.exclude = [0 0 0 0]; 
end
if options.psatz == 1
    options3.psatz=1;
    [sos, Deop] = poslpivar(sos, [nx1 ,nx2],X,d2,options2);
    [sos, De2op] = poslpivar(sos, [nx1 ,nx2],X,d2,options3);
    sos = lpi_eq(sos,Deop+De2op-P); %Dop=Deop+De2op
else
    [sos, Deop] = poslpivar(sos, [nx1 ,nx2],X,d2,options2);
    sos = lpi_eq(sos,Deop-P); %Dop=Deop
end
end
function sos = lpi_ineq_2d(sos,Pop, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sos = lpi_ineq_2d(sos,P, options) sets up inequality constraints such that P>=0
% 
% INPUT
%   prog: SOS program to modify.
%   P: PI dopvar2d variable
%   -Optional INPUTS
%   options: determines if the operator needs to be pure or to include psatz
%   term. 
%   If options.pure(1)=1, then Rxx^o term is excluded. 
%   If options.pure(2)=1, then Ryy^o term is excluded. 
%   If options.pure(3)=1, then R22^oo term is excluded. 
%   If options.pure(4)=1, then R22^ao and R22^bo terms are excluded. 
%   If options.pure(5)=1, then R22^oa and R22^ob terms are excluded. 
%   If options.psatz=1, then the psatz multiplier is included.
% 
% OUTPUT 
%   sos: SOS program
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpi_ineq_2d
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 07_21_2021 
%

% Extract the inputs
switch nargin
    case 1
        error(['Not enough inputs!'])
    case 2
        % If no options are specified, use no special options
        options.psatz = 0;
        options.pure = zeros(1,5);
        options.sep = zeros(1,6);
    case 3
        % If a certain option is not specified, do no use it
        if ~isfield(options,'psatz')
            options.psatz=0;
        end
        if ~isfield(options,'pure')
            options.pure = zeros(1,5);
        end
        if ~isfield(options,'sep')
            options.sep = zeros(1,6);
        end
end

dim = Pop.dim;
if dim(:,1)~=dim(:,2)
    error('Non-symmetric Operators cannot be sign definite. Unable to set the inequality');
end

% We're going to enforce inequality P>=0, by building an operator Deop>=0,
% and setting P=Deop. If psatz is used, build an additional Deop2>=0, and
% set P=Deop+Deop2.

% First, extract the dimensions and domain from P
n0 = dim(1,1);
nx = dim(2,1);
ny = dim(3,1);
n2 = dim(4,1);
dom = Pop.I;

% Next, exclude certain elements from Deop (options2) and Deop2 (options3)
% if specified
options2.exclude = zeros(1,16);
options3.exclude = zeros(1,16);

options2.exclude(2) = options.pure(1);
options3.exclude(2) = options.pure(1); 

options2.exclude(5) = options.pure(2);
options3.exclude(5) = options.pure(2);

tol = 1e-14;
    options3.exclude(8) = 1;
if max(max(Pop.R22{1,1}.C))<tol
    options2.exclude(8) = 1;
else
    options2.exclude(8) = max([options.pure(3),options.pure(4),options.pure(5)]);
    options3.exclude(8) = max([options.pure(3),options.pure(4),options.pure(5)]);
end
if max(max(Pop.R22{2,1}.C))<tol
    options2.exclude(9) = 1;
    options3.exclude(9) = 1;
else
    options2.exclude(9) = options.pure(4);
    options3.exclude(9) = options.pure(4); 
end
if max(max(Pop.R22{3,1}.C))<tol
    options2.exclude(10) = 1;
    options3.exclude(10) = 1;
else
    options2.exclude(10) = options.pure(4);
    options3.exclude(10) = options.pure(4); 
end
if max(max(Pop.R22{1,2}.C))<tol
    options2.exclude(11) = 1;
    options3.exclude(11) = 1;
else
    options2.exclude(11) = options.pure(5);
    options3.exclude(11) = options.pure(5); 
end
if max(max(Pop.R22{1,3}.C))<tol
    options2.exclude(12) = 1;
    options3.exclude(12) = 1;
else
    options2.exclude(12) = options.pure(5);
    options3.exclude(12) = options.pure(5); 
end

% Now, determine appropriate degrees for Deop and Deop2 to match degrees
% for Pop
d = degbalance(Pop);
dx = d.dx;  dy = d.dy;  d2 = d.d2;

% Finally, we can construct the positive operators Deop and Deop2, and use
% them to enforce positivity of P
if options.psatz == 1
    options3.psatz = 1;
    [sos, Deop] = poslpivar_2d(sos, [n0 ,nx, ny, n2],dom,d,options2);
    [sos, De2op] = poslpivar_2d(sos, [n0 ,nx, ny, n2],dom,d,options3);
    sos = lpi_eq_2d(sos,Deop+De2op-Pop); %Pop=Deop+De2op
elseif options.psatz == 2
    options3.psatz = 2;
    [sos, Deop] = poslpivar_2d(sos, [n0 ,nx, ny, n2],dom,d,options2);
    [sos, De2op] = poslpivar_2d(sos, [n0 ,nx, ny, n2],dom,d,options3);
    sos = lpi_eq_2d(sos,Deop+De2op-Pop); %Pop=Deop+De2op
else
    [sos, Deop] = poslpivar_2d(sos, [n0 ,nx, ny, n2],dom,d,options2);
    sos = lpi_eq_2d(sos,Deop-Pop); %Pop=Deop
end

end
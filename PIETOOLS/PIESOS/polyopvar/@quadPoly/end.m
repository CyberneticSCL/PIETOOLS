function e = end(F,k,n)
% A function for indexing quadPoly objects. This will allow operations such
% as obj(end), obj(:,end), obj(1:end,:), etc.
%
% INPUT
% - F:      'quadPoly' object of dimension F.dim = [m,p];
% - k:      integer, the position of the index in which end appears;
% - n:      integer, the number of indices;
%
% OUTPUT
% - e:      integer, m*p for linear indexing (n=1), else m (k=1) or p.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - end
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding SS - 01/16/2026

m = F.dim(1); % find matrix dimensions
p = F.dim(2);

if n == 1  % if obj(end), then use linear indexing
    e = m*p;
else
    if k == 1  % if obj(end,...) then last row
        e = m;
    else       % if obj(..., end) then last column
        e = p;
    end
end
end

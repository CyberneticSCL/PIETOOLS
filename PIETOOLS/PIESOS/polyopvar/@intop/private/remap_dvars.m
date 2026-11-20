function P = remap_dvars(P,m,pos,nT)
% P = REMAP_DVARS(P,M,POS,NT) expresses the sparse coefficient matrix P of
% an M-row matrix of polynomials, with rows (i-1)*(nS+1)+d for
% nS=numel(POS) decision variables (d=1 the constant term, d=l+1 the lth
% decision variable), in terms of a list of NT decision variables, in which
% the lth decision variable of P is decision variable POS(l).
%
% INPUTS
% - P:      m*(nS+1) x c sparse array, for nS=numel(pos), with rows
%           (i-1)*(nS+1)+d holding the coefficients of the polynomial in
%           matrix row i associated with decision variable d;
% - m:      scalar integer specifying the number of matrix rows;
% - pos:    nS x 1 array of integers in {1,...,nT}, with pos(l) the index
%           of the lth decision variable of P in the new list;
% - nT:     scalar integer specifying the number of decision variables in
%           the new list;
%
% OUTPUTS
% - P:      m*(nT+1) x c sparse array, with rows (i-1)*(nT+1)+d holding the
%           coefficients associated with decision variable d of the new
%           list;
%
% NOTES
% - The nonzeros keep their matrix row and column, only the position in the
%   block of rows of the matrix row changes.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - remap_dvars
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
% DJ, 10/15/2026: Initial coding;

nS = numel(pos);
if nS==nT && (nS==0 || isequal(pos(:),(1:nS)'))
    return
end
if m==1 && isequal(pos(:),(1:nS)')
    % Decision variables of P are the first nS: just append empty rows
    P = [P;sparse(nT-nS,size(P,2))];
    return
end
[i,j,v] = find(P);
i = i(:);   j = j(:);   v = v(:);       % row vectors if P has one row
blk = ceil(i/(nS+1));
d = i - (blk-1)*(nS+1);             % 1: constant term, l+1: decision variable l
is_d = d>1;
d_new = ones(size(d));
d_new(is_d) = 1 + pos(d(is_d)-1);
P = sparse((blk-1)*(nT+1)+d_new,j,v,m*(nT+1),size(P,2));

end

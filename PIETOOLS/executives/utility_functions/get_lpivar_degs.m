function Qdeg = get_lpivar_degs(Rop,Top)
% QDEG = GET_LPIVAR_DEGS(ROP,TOP) takes an opvar or dopvar object ROP, and
% returns an array QDEG of degrees to pass to "lpivar" for generating an
% operator variable QOP such that TOP'*QOP can match all monomials
% appearing in ROP, i.e. we can enforce TOP'*QOP==ROP.
%
% INPUT
% - Rop:    'opvar' or 'dopvar' object;
% - Top:    'opvar' or 'dopvar' object;
%
% OUTPUT
% - Qdeg:   1x3 array of nonnegative integers specifying maximal degrees of
%           monomials in call to "lpivar". In particular, calling
%           Qop = lpivar(Top.dim,Qdeg), Qdeg(1) will  determine the degrees
%           of Qop.Q1, Qop.Q2, and Qop.R.R0 in s, Qdeg(2) the degrees of
%           s_dum in Qop.R.R1 and Qop.R.R2, and Qdeg(3) the degrees of s in
%           Qop.R.R2;
%
% NOTES
% The current implementation is very safe but likely also unnecessarily
% expensive. It just checks the maximal degrees of all monomials defining
% Rop, not at all accounting for the operator Top, or distinguishing 
% degrees in s and s_dum. Future updates may improve upon this.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getdegs_lpivar
%
% Copyright (C)2024 PIETOOLS Team
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
% DJ, 01/06/2025: Initial coding;
%

% Check maximal degrees of monomials in Rop.Q1, Rop.Q2, and Rop.R.R0.
Rop_deg1 = 0;
if ~isempty(Rop.Q1) && (isa(Rop.Q1,'polynomial')||isa(Rop.Q1,'dpvar')) && ~isempty(Rop.Q1.degmat)
    Rop_deg1 = max(Rop_deg1,max(Rop.Q1.degmat));
end
if ~isempty(Rop.Q2) && (isa(Rop.Q2,'polynomial')||isa(Rop.Q2,'dpvar')) && ~isempty(Rop.Q2.degmat)
    Rop_deg1 = max(Rop_deg1,max(Rop.Q2.degmat));
end
if ~isempty(Rop.R.R0) && (isa(Rop.R.R0,'polynomial')||isa(Rop.R.R0,'dpvar')) && ~isempty(Rop.R.R0.degmat)
    Rop_deg1 = max(Rop_deg1,max(Rop.R.R0.degmat));
end
% Check maximal degrees of monomials in Rop.R.R1, Rop.R.R2.
Rop_deg2 = 0;
if ~isempty(Rop.R.R1) && (isa(Rop.R.R1,'polynomial')||isa(Rop.R.R1,'dpvar')) && ~isempty(Rop.R.R1.degmat)
    Rop_deg2 = max(Rop_deg2,max(max(Rop.R.R1.degmat)));
end
if ~isempty(Rop.R.R2) && (isa(Rop.R.R2,'polynomial')||isa(Rop.R.R2,'dpvar')) && ~isempty(Rop.R.R2.degmat)
    Rop_deg2 = max(Rop_deg2,max(max(Rop.R.R2.degmat)));
end
Qdeg = [Rop_deg1,Rop_deg2,Rop_deg2-1];  % subtract 1 to account for multiplication with Top

end
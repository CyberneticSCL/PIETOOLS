function [Zth,H] = poly_decomp_1D_econ(Q)
% poly_decomp_1D_econ creates a monomial-matrix decomposition of the matrix-valued polynomial
% Q in a single spatial variable such that Q(s)=H*Zth(s) where Z is the
% smallest vector of n-dimensional monomial bases needed
%
% Note: This representation is unique for the monomial basis Zth. The use 
% of alternative bases for the polynomial matrices is not currently
% supported.
%
% Note: To get a left decomposition, use [ZthT,HT] = poly_decomp_1D_econ(Q') and 
% then Q=ZthT'*HT'
%
% INPUTS 
%   Q: A pvar object of dimension n times m with a single independent variable (there is currently no error
%   checking to ensure the absence of other independent variables.
%
% OUTPUT 
%   Zth - a pvar object of dimension q times m which is a matrix of vector-valued monomial bases
%   H   - a non-sparse matrix of dimension n times q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poly_decomp_1D_econ
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
% NOTES:
% Distributed with DelayTOOLS
% Compatable with MULTIPOLY and SOSTOOLS as of June 2013
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial Coding:
% 6/26/20 - MMP  -  titled poly_decomp_1D_econ
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/20/2026: Column monomials from the degmat (no one-argument monomials
%                   on the path); a polynomial without variables as a constant.

% A polynomial without variables (multipoly drops s1 from A+0*s1) failed    % MMP, 11/20/2026
% at Q.Varname(1) below; it is a constant, with H as 'double'               % MMP, 11/20/2026
% if ~isa(Q,'polynomial')                                                   % MMP, 11/20/2026 (was)
if ~isa(Q,'polynomial') || isempty(Q.varname)                               % MMP, 11/20/2026
    [nd md]=size(Q);
    Zth=polynomial(eye(md));
%   H = Q;                                                                  % MMP, 11/20/2026 (was)
    H = double(Q);                                                          % MMP, 11/20/2026
    return
end
if length(Q.Varname)>1
    error('Q has more than 1 independent variable')
end
ss = Q.Varname(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first step is to construct the separable representation of the
% Polynomial Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Q, we need to first convert to the forms
% Q=H Z(s) or for us this is Q=HZ(theta)
% however, we need to be consistent about our Z

Qt=Q;
[nd md]=size(Qt);
bZth=[];
for i=1:md
%   sZth{i}=monomials(Qt(:,i));                                             % MMP, 11/20/2026 (was)
    % Monomials of column i, sorted by degree, from its degmat, as the      % MMP, 11/20/2026
    % one-argument monomials(p) of multipoly did; the monomials on the path % MMP, 11/20/2026
    % (SOSTOOLS) needs variables and degrees ("Not enough input arguments") % MMP, 11/20/2026
    Qi = Qt(:,i);   Nt = size(Qi.degmat,1);                                 % MMP, 11/20/2026
    sZth{i}=polynomial(speye(Nt),sortrows(Qi.degmat),Qi.varname,[Nt 1]);    % MMP, 11/20/2026
    bZth=blkdiag(bZth,sZth{i});
    nZ(i)=length(sZth{i});
end


% Now use bigZth to find H!
for i=1:nd
    for j=1:md
        [CQij,Ztemp,etemp] = poly2basis(Qt(i,j),sZth{j}); % examine each element of NQ
        bigC(i,(sum(nZ(1:(j-1)))+1):sum(nZ(1:j)))=CQij.';
    end
end
H=full(bigC);
Zth=bZth;

% Q'=H*Zth so Q=Zth'*H'
% Zth*H-Q %uncomment to verify the representation (should be 0)

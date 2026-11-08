function [Aout,Bout] = lr_multiply(L,A,B,R)                                 % MMP, 09/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [AOUT,BOUT] = LR_MULTIPLY(L,A,B,R) returns, for the matrix C(d) with
% vec(C(d)) = A + B'*d, the coefficients of L*C(d)*R, using
% vec(L*C*R) = kron(R.',L)*vec(C):
%
%       vec(L*C(d)*R) = Aout + Bout'*d,
%       Aout = vec(L*reshape(A,size(L,2),size(R,1))*R),
%       Bout = B*kron(R.',L).'.
%
% INPUT
% - L:      'double', p x mr, sparse or full;
% - A:      'double', (mr*nc) x 1, the vec of the constant term of C(d);
% - B:      'double', nd x (mr*nc), row i the vec of the coefficient of
%           d_i, as a row; not used, and may be [], when BOUT is not
%           requested;
% - R:      'double', nc x r, sparse or full.
%
% OUTPUT
% - Aout:   'double', (p*r) x 1;
% - Bout:   'double', nd x (p*r).
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lr_multiply
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
% MMP, 09/26/2026: Optional fifth input K, the product kron(R.',L).' formed
%                  once by the caller; B skipped when only Aout is requested
%                  ('unpack_sheets' maps B itself).
% MMP, 09/30/2026: K removed: no caller passed it any more. Bout is
%                  B*kron(R.',L).' again, the product the K default formed,
%                  reached on the unfused path of 'unpack_sheets' only.
% MMP, 10/28/2026: Header to the header protocol, step comments in the body;
%                  code unchanged. Deleted: the cost rationale for forming
%                  the Kronecker product and the help of the input K. No
%                  initial-coding entry: MP, 08/22/2026, as a local function
%                  of 'possopvar' (since renamed 'sopquadvar'), moved to
%                  private/ by MMP, 09/21/2026. Earlier signatures:
% function [Aout,Bout] = lr_multiply(L,A,B,R,K)                             % MMP, 09/30/2026 (was)
% function [Aout,Bout] = lr_multiply(L,A,B,R)                               % MMP, 09/26/2026 (was)
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed;
%                  code unchanged.

% Constant term: A as the mr x nc matrix, multiplied, vectorized
X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;
Aout = Y(:);
if nargout<2,  return;  end     % B unused: no Kronecker product            % MMP, 09/26/2026

% Decision part: every row of B at once, times kron(R.',L).'
% Bout = B*kron(R.',L).';                                                   % MMP, 09/26/2026 (was)
% if nargin<5 || isempty(K),  K = kron(R.',L).';  end                       % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% Bout = B*K;                                                               % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
Bout = B*kron(R.',L).';     % the product the K default formed              % MMP, 09/30/2026

end

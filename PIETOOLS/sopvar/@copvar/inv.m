function [Pinv,info] = inv(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PINV,INFO] = INV(P[,OPTS]) the inverse of a square 1-D 'copvar' over the
% spaces R^k and L2^m[s] (either may be absent), the container counterpart
% of @opvar/inv. With the grid written as
%   P = [Pm  Q1;  Q2  R],   Pm: R^k -> R^k,  Q1: L2 -> R^k,  Q2: R^k -> L2,
% R: L2 -> L2, the L2 block is inverted by @sopvar/inv (Gohberg-Krein) and
% the rest is the block inverse through the finite-dimensional Schur
% complement (Shivakumar, Das, Peet, arXiv 2208.13104, Lemma 18):
%   T = Pm - Q1 Rh Q2  (k x k),   Rh = R^{-1},
%   P^{-1} = [T^{-1},        -T^{-1} Q1 Rh;
%             -Rh Q2 T^{-1},  Rh + Rh Q2 T^{-1} Q1 Rh].
% Only one operator inverse (Rh) and one matrix inverse are taken; the stock
% 'inv_opvar' takes the other Schur complement, (R - Q2 Pm^{-1} Q1)^{-1},
% which is a second operator inverse of a modified kernel. The products are
% block compositions of the class ('sopvar' mtimes), exact on polynomials.
%
% INPUTS
% - P:      'copvar', square: the same spaces and dimensions in and out, at
%           most one R^k and one L2[s] space, one spatial variable;
% - opts:   (optional) the options of @sopvar/inv (N, tol, deg0, degmax,
%           deg, tol_rank).
% OUTPUTS
% - Pinv:   'copvar' over the same spaces in the same order as P;
% - info:   the INFO of @sopvar/inv (empty when there is no L2 space), plus
%           condT, the condition number of the Schur complement, when both
%           spaces are present.
%
% NOTES
% A structurally zero block ([]) of P is respected: a zero Q1 or Q2 leaves
% T = Pm and the off-diagonal blocks of the inverse zero. The R^k -> R^k
% block Pm may itself be [], in which case T = -Q1 Rh Q2.
% N-D containers are refused: the Gohberg-Krein construction is one
% dimensional, and no inverse formula for a non-separable 2-D kernel exists
% in the toolbox ('opvar2d' falls back to least squares, 'mrdivide').
%
% See also SOPVAR/INV, GETCONTROLLER_SOP, GETOBSERVER_SOP, COPVAR2OPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - inv (copvar)
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
% Initial coding MMP, 10/09/2026.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,    opts = struct();    end
if ~isa(P,'copvar'),    error('copvar:inv:class','Input must be a copvar.');   end
if numel(P.vars)>1
    error('copvar:inv:nd','INV is implemented for one spatial variable; this container has %d (%s).',numel(P.vars),strjoin(P.vars,','));
end
[M,Nc] = size(P);
if M~=Nc || ~isequal(P.space_out,P.space_in) || ~isequal(P.dim_out(:),P.dim_in(:))
    error('copvar:inv:square','The container must map a space list to itself (same spaces and dimensions in and out).');
end
isL2 = any(P.space_out,2);
l = find(isL2);     r = find(~isL2);
if numel(l)>1 || numel(r)>1
    error('copvar:inv:spaces','At most one R^k and one L2[s] space are supported; the grid is %d x %d.',M,Nc);
end
C = cell(M,M);
info = struct();
if isempty(l)
    % R^k -> R^k alone
    k = P.dim_out(r);
    C{r,r} = matblock(blockmat(P.C{r,r},k)\eye(k),k);
elseif isempty(r)
    [C{l,l},info] = inv(P.C{l,l},opts);
else
    k = P.dim_out(r);
    [Rh,info] = inv(P.C{l,l},opts);
    Q1 = P.C{r,l};  Q2 = P.C{l,r};
    T = blockmat(P.C{r,r},k);
    if ~isempty(Q1) && ~isempty(Q2)
        T = T - blockmat(Q1*Rh*Q2,k);
    end
    info.condT = cond(T);
    Th = matblock(T\eye(k),k);
    C{r,r} = Th;
    if ~isempty(Q1),    C{r,l} = -(Th*(Q1*Rh));     end
    if ~isempty(Q2),    C{l,r} = -((Rh*Q2)*Th);     end
    if ~isempty(Q1) && ~isempty(Q2)
        C{l,l} = Rh + (Rh*Q2)*(Th*(Q1*Rh));
    else
        C{l,l} = Rh;
    end
end
Pinv = copvar(C);
end


function Mat = blockmat(B,k)
% the k x k matrix of an R^k -> R^k block, [] being zero
if isempty(B),  Mat = zeros(k);     return,     end
Mat = full(B.params{1});
if ~isequal(size(Mat),[k k])
    error('copvar:inv:matblock','An R^%d -> R^%d block stores a %d x %d parameter.',k,k,size(Mat,1),size(Mat,2));
end
end

function B = matblock(Mat,k)
% the 'sopvar' block of a k x k matrix on R^k, as 'opvar2sopvar' builds it
vars = struct;  vars.in = {};   vars.out = {};
dom = struct;   dom.in = zeros(0,2);    dom.out = zeros(0,2);
B = sopvar({Mat},vars,cell(1,0),cell(1,0),dom,[k k]);
end

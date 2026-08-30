function y = apply_nopvar_quad(P,xfun,s,nq)
% Y = APPLY_NOPVAR_QUAD(P,XFUN,S,NQ) evaluates (P*x)(s) directly from the
% definition of the 'nopvar' class, using Gauss--Legendre quadrature. For
% polynomial x the result is exact provided NQ is large enough.
%
% INPUTS
% - P:      'nopvar' object with N spatial directions and dimension [m,n];
% - xfun:   function handle mapping an Nx1 point t to an n x 1 vector;
% - s:      Nx1 array specifying the point at which to evaluate (P*x);
% - nq:     scalar integer specifying the number of quadrature nodes to use
%           along each integrated direction (default 20);
%
% OUTPUTS
% - y:      m x 1 array specifying the value of (P*x) at s;
%
% NOTES
% This routine deliberately implements nothing more than the definition
%   (P*x)(s) = sum_{j} int I_j(s,t) (Im o Zd(s))^T Cj (In o Zj(t)) x(t) dt,
% so that it can serve as an independent reference for the (much faster)
% coefficient-level operations implemented in the class.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - apply_nopvar_quad
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
% MMP, 08/23/2026: Initial coding

if nargin<4
    nq = 20;
end
deg = P.deg(:);
N   = numel(deg);
dom = P.dom;
dm  = P.dim;
m = dm(1);      n = dm(2);

% Monomial vector Zd(s) = Zd1(s1) o ... o ZdN(sN)
Zs = 1;
for i=1:N
    Zs = kron(Zs,s(i).^(0:deg(i))');
end
Lfac = kron(speye(m),Zs');

y = zeros(m,1);
sz_C = [size(P.C),1];   sz_C = [sz_C(1:N),1];
idcs = cell(1,N);
sub  = cell(1,N);
for ii=1:numel(P.C)
    C = P.C{ii};
    if isempty(C) || nnz(C)==0
        continue
    end
    [idcs{:}] = ind2sub(sz_C,ii);
    jj = cell2mat(idcs)-1;                  % 0, 1 or 2 along each direction

    % Quadrature nodes and weights along each direction: a multiplier
    % corresponds to evaluation at t=s, rather than integration
    nod = cell(N,1);    wgt = cell(N,1);
    for i=1:N
        if jj(i)==0
            nod{i} = s(i);      wgt{i} = 1;
        elseif jj(i)==1
            [nod{i},wgt{i}] = gauss_leg(nq,dom(i,1),s(i));
        else
            [nod{i},wgt{i}] = gauss_leg(nq,s(i),dom(i,2));
        end
    end

    npts = cellfun(@numel,nod);
    for kk=1:prod(npts)
        [sub{:}] = ind2sub([npts(:)',1],kk);
        idx = cell2mat(sub);
        t = zeros(N,1);     w = 1;
        for i=1:N
            t(i) = nod{i}(idx(i));
            w = w*wgt{i}(idx(i));
        end
        % Only the integrated directions contribute dummy monomials
        Zt = 1;
        for i=1:N
            if jj(i)>0
                Zt = kron(Zt,t(i).^(0:deg(i))');
            end
        end
        y = y + w*((Lfac*C*kron(speye(n),Zt))*xfun(t));
    end
end

end



%% Gauss--Legendre nodes and weights on [a,b], via the Golub--Welsch method
function [x,w] = gauss_leg(nn,a,b)

if b==a
    x = a;      w = 0;
    return
end
bet = 0.5./sqrt(1-(2*(1:nn-1)).^(-2));
T = diag(bet,1)+diag(bet,-1);
[V,D] = eig(T);
[x,ord] = sort(diag(D));
w = 2*V(1,ord).^2;
x = (a+b)/2 + (b-a)/2*x;
w = (b-a)/2*w(:);

end

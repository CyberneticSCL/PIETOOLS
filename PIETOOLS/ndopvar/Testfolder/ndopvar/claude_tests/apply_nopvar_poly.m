function y = apply_nopvar_poly(P,Xc,dx,s)
% Y = APPLY_NOPVAR_POLY(P,XC,DX,S) evaluates (P*x)(s) exactly from the
% definition of the 'nopvar' class, for a polynomial input x, performing all
% integrals in closed form.
%
% INPUTS
% - P:      'nopvar' object with N spatial directions and dimension [m,n];
% - Xc:     K x n array of coefficients defining
%               x_k(t) = sum_{alpha} Xc(alpha,k)*t^alpha,
%           with the exponents alpha enumerated as in EXPGRID(DX), i.e. in
%           the same Kronecker ordering used by the class;
% - dx:     Nx1 array specifying the componentwise degree of x;
% - s:      Nx1 array specifying the point at which to evaluate (P*x);
%
% OUTPUTS
% - y:      m x 1 array specifying the value of (P*x) at s;
%
% NOTES
% As with APPLY_NOPVAR_QUAD, this routine implements only the definition of
% the operator, and is intended as an independent reference for testing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - apply_nopvar_poly
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

deg = P.deg(:);
N   = numel(deg);
dom = P.dom;
dm  = P.dim;
m = dm(1);      n = dm(2);
dx = dx(:);

Ex = expgrid(dx);               % K x N array of exponents appearing in x
K  = size(Ex,1);

Zs = 1;
for i=1:N
    Zs = kron(Zs,s(i).^(0:deg(i))');
end
Lfac = kron(speye(m),Zs');

y = zeros(m,1);
sz_C = [size(P.C),1];   sz_C = [sz_C(1:N),1];
idcs = cell(1,N);
for ii=1:numel(P.C)
    C = P.C{ii};
    if isempty(C) || nnz(C)==0
        continue
    end
    [idcs{:}] = ind2sub(sz_C,ii);
    jj = cell2mat(idcs)-1;                  % 0, 1 or 2 along each direction
    int_dims = find(jj>0);
    Eb = expgrid(deg(int_dims));            % exponents of the dummy basis
    Kb = size(Eb,1);

    % W(beta,alpha) = prod_{i: multiplier} s_i^{alpha_i}
    %                   * prod_{i: integral} int_{D_i(s)} t^{beta_i+alpha_i}
    W = ones(Kb,K);
    for i=1:N
        if jj(i)==0
            W = W.*repmat((s(i).^Ex(:,i))',[Kb,1]);
        else
            g = Eb(:,int_dims==i) + Ex(:,i)';
            if jj(i)==1
                W = W.*((s(i).^(g+1) - dom(i,1).^(g+1))./(g+1));
            else
                W = W.*((dom(i,2).^(g+1) - s(i).^(g+1))./(g+1));
            end
        end
    end
    % Block k of the integrated vector (In o Zj(t))*x(t) is Zj(t)*x_k(t)
    y = y + Lfac*C*reshape(W*Xc,[Kb*n,1]);
end

end



%% Enumerate the exponents of Zd1(s1) o ... o ZdN(sN)
function E = expgrid(degvec)

degvec = degvec(:);
E = zeros(1,0);
for i=1:numel(degvec)
    z = (0:degvec(i))';
    E = [repelem(E,numel(z),1), repmat(z,size(E,1),1)];
end

end

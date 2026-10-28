function [Tx,Fop,Kop,info] = compute_PIETOOLS_PDE_periodic_map(PDE,Pop_f2x,Pop_b2x,Gop,Hop)
% Reconstruct a 1D PDE with singular homogeneous BCs and mixed spatial orders.
% The existing FTC expansion gives x = Pop_b2x*xb + Pop_f2x*xf and BCs
% G*xb + Hop*xf = 0, where xf takes the declared derivative of each state.
% Retain additional moment states v0=Fop*x, and return
% x=Tx*[v0;xf], Kop*[v0;xf]=0.
% Implements Definition 1, Lemma 5, Proposition 6 and Lemma 7 of
% Jagt, Chernyshenko and Peet, arXiv:2503.22896v5. Hop uses the negative
% of the paper's H convention. Q preserves the converter's core-boundary
% ordering, including separate vector-valued PDE components.
%
% Initial coding: DB 10/27/2026.

s = PDE.vars(1,1); theta = PDE.vars(1,2);
a = PDE.dom(1,1); b = PDE.dom(1,2); L = b-a;
G = double(Gop.P); Q = polynomial(Pop_b2x.Q2);
H = polynomial(Hop.Q1);
n = sum(PDE.x_tab(:,2)); nb = size(G,2);

%% Lemma 3 
J = eye(nb);
tol = eps;

%% Gauss-Jordan elimination with partial row pivoting. 
pivots = zeros(1,0); r = 0;
for j=1:nb
    % Search for the pivot, the row with largest entry in column j
    [pivot,idx] = max(abs(G(r+1:end,j)));
    % If this column has no significant entry in the remaining rows,
    % it cannot supply another independent pivot.
    if isempty(pivot) || pivot<=tol
        continue
    end
    % The largest entry of column j is at row idx. Sweep row r with idx
    idx = idx+r; 
    r = r+1;
    G([r,idx],:) = G([idx,r],:); 
    J([r,idx],:) = J([idx,r],:);
   % Normalize the pivot row 'r' so its entry in column j is 1.
    pivot = G(r,j);
    G(r,:) = G(r,:)/pivot; 
    J(r,:) = J(r,:)/pivot;
    % subtract to eliminate column j from the other i rows (i \neq r) by the r-th row times G(i,j)
    % Setting factors(r)=0 leaves the normalized pivot row unchanged.
    factors = G(:,j); factors(r) = 0;
    G = G-factors*G(r,:); 
    J = J-factors*J(r,:);
   % Record the pivot column and stop when every row has a pivot.
    pivots(end+1) = j;
    if r==nb, 
        break; 
    end
end
% Nonpivot columns identify the core-boundary states.
%  m = nb-r is the nullity of G and determines how many
% auxiliary average states are needed.
%returns the values in {1,...,nb} which are not in pivots, in the order that they
% appear in {1,...,nb}
free = setdiff(1:nb,pivots,'stable');
m = numel(free);
% This routine is called for the rank-deficient case. If m=0, something is
% going wrong.
if m==0
    error('Boundary rank is numerically ambiguous. Check the boundary conditions.')
end
G1 = G(1:r,:); 
H = J*H;

%% For spatial derivatives higher then second order
% Each free FTC column is a single scalar state's Taylor monomial. Record
% its component and degree, including traces above the first derivative.
source_components = zeros(m,1);
source_scalars = zeros(m,1);
trace_orders = zeros(m,1);
ends = cumsum(PDE.x_tab(:,2));
for j=1:m
    Qcol = Q(:,free(j));
    source_scalars(j) = find(any(Qcol.coefficient,1),1);
    source_components(j) = find(source_scalars(j)<=ends,1);
    trace_orders(j) = max(full(sum(Qcol.degmat,2)));
end

%% Lemma 5
% Preserve the second-order average choice. For other orders, try averages
% from the free Taylor columns of the existing FTC basis.
is_second_order = all(PDE.x_tab(:,4)==2);
if is_second_order
    V = double(subs(Q,s,a)); 
    D = double(diff(Q,s));
    F3 = V(:,free)'/L + 2*(s-a)*D(:,free)'/L^2;
else
    F3 = Q(:,free)'/L;
end
G3 = double(int(F3*Q,s,a,b));
if rank([G1;G3])<nb
    if is_second_order
        F3 = 2*(2*b+a-3*s)*V(:,free)'/L^2 ...
            -6*(b+a-2*s)*D(:,free)'/L^3;
    else
        % Use polynomial moments dual to the free endpoint traces:
        % int_a^b F3(s)*Q(s) ds selects exactly the free columns.
        % A shifted Legendre expansion avoids solving a Taylor Gram matrix.
        max_order = max(PDE.x_tab(:,4));
        P = cell(max_order,1);
        P{1} = polynomial(1);
        z = (2*s-a-b)/L;
        if max_order>1
            P{2} = z; 
        end
        for ell=1:max_order-2
            P{ell+2} = ((2*ell+1)*z*P{ell+1}-ell*P{ell})/(ell+1);
        end
        F3 = polynomial(zeros(m,n));
        for j=1:m
            k = trace_orders(j);
            d = PDE.x_tab(source_components(j),4);
            for ell=k:d-1
                weight = (2*ell+1)*(-1)^(ell+k)*prod(ell-k+1:ell+k) ...
                    /(factorial(k)*L^(k+1));
                F3(j,source_scalars(j)) = F3(j,source_scalars(j)) ...
                    + weight*P{ell+1};
            end
        end
    end
    G3 = double(int(F3*Q,s,a,b));
end
Ghat = [G1;G3];
if rank(Ghat)<nb
    error('Unable to construct independent average states for the singular boundary conditions.')
end
%% Compose the moment functional with the already constructed FTC map.
Fop = opvar(); Fop.I = PDE.dom;
Fop.var1 = s; Fop.var2 = theta; Fop.dim = [m,0;0,n];
Fop.Q1 = F3;
Fftc = Fop*Pop_f2x;
H3 = subs(Fftc.Q1,s,theta);
Ginv = Ghat\eye(nb);
C = -Q*Ginv*[subs(H(1:r,:),s,theta); H3];
T0 = Q*Ginv*[zeros(r,m);eye(m)];

Tx = opvar(); Tx.I = PDE.dom;
Tx.var1 = s; Tx.var2 = theta; Tx.dim = [0,m;n,n];
Tx.Q2 = T0;
Tx.R.R0 = Pop_f2x.R.R0;
Tx.R.R1 = Pop_f2x.R.R1+C; Tx.R.R2 = Pop_f2x.R.R2+C;

Kop = opvar(); Kop.I = PDE.dom;
Kop.var1 = s; Kop.var2 = theta; Kop.dim = [m,m;0,n];
Kop.Q1 = H(r+1:end,:);

%% Save conversion info
info = struct('reconstruction',Tx,'n_aux',m, ...
    'boundary_rank',r,'jordan_transform',J,'free_boundary_columns',free, ...
    'source_components',source_components);
end

function M = mixed_domain_colloc(name,p,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M = MIXED_DOMAIN_COLLOC(NAME,P,N) discretizes the PDE of a system of
% MIXED_DOMAIN_SYSTEMS by Chebyshev collocation on the PDE itself, as a
% reference independent of the parser, the PIE and PIESIM. Every state is
% a vector of values at the Chebyshev-Lobatto nodes of [0,1] (N+1 per
% direction, s1 fastest in 2-D). Each equation holds at the nodes off the
% boundary, each BC at the nodes of its edge (a corner node takes the BC of
% its s1 edge). With u = [u_I; u_B] (dynamic and boundary values), the BC
% rows give u_B = G u_I and the equations u_I' = L u_I.
% M = MIXED_DOMAIN_COLLOC('selftest',[],N) checks the nodes, derivative and
% quadrature on x^3.
%
% INPUT
% - name:   char, one of 'S2_ode_heat2d_edge', 'S3_edge1d_heat2d',
%           'S4_wentzell', 'S5_chain_ode_1d_2d', 'S6_indomain_1d_2d',
%           'S7_ode_2d_indomain', 'S8_two_edges', 'S9_ode_2d_flux',
%           'C1_heat2d_scalar_w', or 'selftest';
% - p:      (optional) 'struct' of parameters, as in MIXED_DOMAIN_SYSTEMS;
% - N:      (optional) double, the polynomial degree per direction
%           (default 24).
%
% OUTPUT
% - M:      'struct' with fields
%             L, G       double, the reduced generator and boundary map;
%             idxI, idxB double column, the global indices of u_I and u_B;
%             nodes      (N+1) x 1 double, the nodes, ascending;
%             w          (N+1) x 1 double, the Clenshaw-Curtis weights;
%             state      struct array, per state: kind ('ode', 's1', 's2'
%                        or 's12'), idx (global indices of its values) and
%                        shape;
%             eig        eigenvalues of L, real part descending;
%             condBB     condition number of the boundary block;
%           for 'selftest', fields eD, eD2, eW (errors of the derivative,
%           second derivative and quadrature) and asc (true if the nodes
%           ascend from 0 to 1).
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mixed_domain_colloc
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
% Initial coding MMP, 11/12/2026

if nargin<2 || isempty(p), p = struct(); end
if nargin<3, N = 24; end
% Nodes, differentiation matrices and quadrature weights on [0,1]
[x,Dm] = cheb01(N);  D2m = Dm*Dm;  w = ccw(N);
if strcmp(name,'selftest')
    M = struct('eD',max(abs(Dm*x.^3-3*x.^2)),'eD2',max(abs(D2m*x.^3-6*x)), ...
               'eW',abs(w'*x.^3-1/4),'asc',all(diff(x)>0)&&x(1)==0&&abs(x(end)-1)<1e-15);
    return
end
n = N+1;  I = speye(n);
e0 = 1;  e1 = n;                           % node indices of s=0 and s=1
% The states of each system and their global indices
switch name
    case {'S3_edge1d_heat2d','S4_wentzell','S6_indomain_1d_2d'}
        kinds = {'s1','s12'};
    case 'S5_chain_ode_1d_2d'
        kinds = {'ode','s1','s12'};
    case 'S8_two_edges'
        kinds = {'s1','s2','s12'};
    case {'S9_ode_2d_flux','S2_ode_heat2d_edge','S7_ode_2d_indomain'}
        kinds = {'ode','s12'};
    case 'C1_heat2d_scalar_w'
        kinds = {'s12'};
    otherwise
        error('mixed_domain_colloc:name','No collocation model for %s.',name);
end
st = struct('kind',{},'idx',{},'shape',{});
nglob = 0;
for k=1:numel(kinds)
    switch kinds{k}
        case 'ode',       sh = [1 1];
        case {'s1','s2'}, sh = [n 1];
        case 's12',       sh = [n n];
    end
    st(k) = struct('kind',kinds{k},'idx',nglob+(1:prod(sh)),'shape',sh);
    nglob = nglob+prod(sh);
end
% Full operator: the equation in every row, then BC rows replace the rows of
% boundary nodes (isB). 2-D operators act on u(:) with s1 fastest.
Lf = sparse(nglob,nglob);  isB = false(nglob,1);
D1_2 = kron(I,Dm);  D11_2 = kron(I,D2m);  D2_2 = kron(Dm,I);  D22_2 = kron(D2m,I);
lin2 = @(i1,i2) (i2-1)*n+i1;
all1 = (1:n)';
edge = @(which) edge_nodes(which,n,lin2);
getp = @(f,d) getpar(p,f,d);
switch name
    case 'S3_edge1d_heat2d'
        r = getp('r',0);
        a = st(1).idx;  b = st(2).idx;
        Lf(a,a) = D2m+r*I;  Lf(b,b) = D11_2+D22_2;
        [Lf,isB] = bc_dir1(Lf,isB,a(e0));  [Lf,isB] = bc_dir1(Lf,isB,a(e1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=0'),[]);  [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=1'),[]);
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=0 nc'),a(2:n-1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=1 nc'),[]);
    case 'S4_wentzell'
        c0 = getp('c',1);  r = getp('r',0);
        a = st(1).idx;  b = st(2).idx;
        Lf(a,a) = D2m+r*I;
        Lf(a,b) = c0*D2_2(lin2(all1,1),:);         % c x2_s2(s1,0) at each s1 node
        Lf(b,b) = D11_2+D22_2+r*speye(n^2);
        [Lf,isB] = bc_dir1(Lf,isB,a(e0));  [Lf,isB] = bc_dir1(Lf,isB,a(e1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=0'),[]);  [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=1'),[]);
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=0 nc'),a(2:n-1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=1 nc'),[]);
    case 'S5_chain_ode_1d_2d'
        a0 = getp('a',1);
        X = st(1).idx;  a = st(2).idx;  b = st(3).idx;
        Lf(X,X) = -a0;  Lf(a,a) = D2m;  Lf(b,b) = D11_2+D22_2;
        % x1'(0) = X, x1(1) = 0
        Lf(a(e0),:) = 0;  Lf(a(e0),a) = Dm(e0,:);  Lf(a(e0),X) = -1;  isB(a(e0)) = true;
        [Lf,isB] = bc_dir1(Lf,isB,a(e1));
        % x2_s1(0,s2) = (1-s2) X on every node of the edge, corners included
        en = edge('s1=0');
        for j=1:numel(en)
            row = b(en(j));  [~,i2] = ind2sub([n n],en(j));
            Lf(row,:) = 0;  Lf(row,b) = D1_2(en(j),:);  Lf(row,X) = -(1-x(i2));  isB(row) = true;
        end
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=1'),[]);
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=0 nc'),a(2:n-1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=1 nc'),[]);
    case 'S6_indomain_1d_2d'
        g = getp('g',1);
        a = st(1).idx;  b = st(2).idx;
        Lf(a,a) = D2m;  Lf(a,b) = g*kron(w(:)',I);  % g int x2 ds2 at each s1 node
        Lf(b,b) = D11_2+D22_2;  Lf(b,a) = g*kron(ones(n,1),I);
        [Lf,isB] = bc_dir1(Lf,isB,a(e0));  [Lf,isB] = bc_dir1(Lf,isB,a(e1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=0'),[]);  [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=1'),[]);
        [Lf,isB] = bc_neu2(Lf,isB,b,edge('s2=0 nc'),D2_2);  [Lf,isB] = bc_neu2(Lf,isB,b,edge('s2=1 nc'),D2_2);
    case 'S8_two_edges'
        r3 = getp('r3',0);
        a = st(1).idx;  c = st(2).idx;  b = st(3).idx;
        Lf(a,a) = D2m;  Lf(c,c) = D2m+r3*I;  Lf(b,b) = D11_2+D22_2;
        [Lf,isB] = bc_dir1(Lf,isB,a(e0));  [Lf,isB] = bc_dir1(Lf,isB,a(e1));
        [Lf,isB] = bc_dir1(Lf,isB,c(e0));  [Lf,isB] = bc_dir1(Lf,isB,c(e1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=0'),[]);
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=1'),c(:));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=0 nc'),a(2:n-1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=1 nc'),[]);
    case 'S9_ode_2d_flux'
        a0 = getp('a',0.5);  c0 = getp('c',1);
        X = st(1).idx;  b = st(2).idx;
        Lf(X,X) = a0;  Lf(X,b) = c0*w(:)'*D2_2(lin2(all1,1),:);   % c int x2_s2(s1,0) ds1
        Lf(b,b) = D11_2+D22_2;
        [Lf,isB] = bc_neu2(Lf,isB,b,edge('s1=0'),D1_2);  [Lf,isB] = bc_neu2(Lf,isB,b,edge('s1=1'),D1_2);
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=0 nc'),repmat(X,n-2,1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=1 nc'),[]);
    case 'S2_ode_heat2d_edge'
        a0 = getp('a',1);
        X = st(1).idx;  b = st(2).idx;
        Lf(X,X) = -a0;  Lf(b,b) = D11_2+D22_2;
        [Lf,isB] = bc_neu2(Lf,isB,b,edge('s1=0'),D1_2);  [Lf,isB] = bc_neu2(Lf,isB,b,edge('s1=1'),D1_2);
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=0 nc'),repmat(X,n-2,1));
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=1 nc'),[]);
    case 'S7_ode_2d_indomain'
        a0 = getp('a',1);  bb = getp('b',1);  g = getp('g',0.5);
        X = st(1).idx;  b = st(2).idx;
        Lf(X,X) = -a0;  Lf(X,b) = g*kron(w(:)',w(:)');
        Lf(b,b) = D11_2+D22_2-bb*speye(n^2);  Lf(b,X) = g;
        [Lf,isB] = bc_neu2(Lf,isB,b,edge('s1=0'),D1_2);  [Lf,isB] = bc_neu2(Lf,isB,b,edge('s1=1'),D1_2);
        [Lf,isB] = bc_neu2(Lf,isB,b,edge('s2=0 nc'),D2_2);  [Lf,isB] = bc_neu2(Lf,isB,b,edge('s2=1 nc'),D2_2);
    case 'C1_heat2d_scalar_w'
        b = st(1).idx;
        Lf(b,b) = D11_2+D22_2;
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=0'),[]);  [Lf,isB] = bc_dir2(Lf,isB,b,edge('s1=1'),[]);
        [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=0 nc'),[]);  [Lf,isB] = bc_dir2(Lf,isB,b,edge('s2=1 nc'),[]);
end
% Eliminate the boundary values: L_BB u_B + L_BI u_I = 0
iB = find(isB);  iI = find(~isB);
LBB = full(Lf(iB,iB));  LBI = full(Lf(iB,iI));
G = -(LBB\LBI);
L = full(Lf(iI,iI))+full(Lf(iI,iB))*G;
M = struct('L',L,'G',G,'idxI',iI,'idxB',iB,'nodes',x,'w',w,'state',st,'N',N);
ev = eig(L);
[~,o] = sort(real(ev),'descend');
M.eig = ev(o);
M.condBB = cond(LBB);
end



function [Lf,isB] = bc_dir1(Lf,isB,row)
% [LF,ISB] = BC_DIR1(LF,ISB,ROW) sets the row ROW to u(ROW) = 0.
Lf(row,:) = 0;  Lf(row,row) = 1;  isB(row) = true;
end

function [Lf,isB] = bc_dir2(Lf,isB,b,en,src)
% [LF,ISB] = BC_DIR2(LF,ISB,B,EN,SRC) sets u = src on the edge nodes EN of
% the 2-D state with global indices B; SRC holds the global index of the
% source value at each edge node, or is [] for u = 0.
for j=1:numel(en)
    row = b(en(j));
    Lf(row,:) = 0;  Lf(row,row) = 1;
    if ~isempty(src), Lf(row,src(j)) = -1; end
    isB(row) = true;
end
end

function [Lf,isB] = bc_neu2(Lf,isB,b,en,Dk)
% [LF,ISB] = BC_NEU2(LF,ISB,B,EN,DK) sets (DK u)(node) = 0 on the edge
% nodes EN of the 2-D state with global indices B.
for j=1:numel(en)
    row = b(en(j));
    Lf(row,:) = 0;  Lf(row,b) = Dk(en(j),:);  isB(row) = true;
end
end

function en = edge_nodes(which,n,lin2)
% EN = EDGE_NODES(WHICH,N,LIN2) returns the linear indices of the nodes of
% an edge of the n x n grid; 'nc' excludes the two corners.
all1 = (1:n)';  inn = (2:n-1)';
switch which
    case 's1=0',     en = lin2(1,all1);
    case 's1=1',     en = lin2(n,all1);
    case 's2=0 nc',  en = lin2(inn,1);
    case 's2=1 nc',  en = lin2(inn,n);
end
en = en(:);
end

function v = getpar(p,f,d)
% V = GETPAR(P,F,D) returns the field F of the 'struct' P, or D if absent.
if isfield(p,f), v = p.(f); else, v = d; end
end

function [x,D] = cheb01(N)
% [X,D] = CHEB01(N) returns the N+1 Chebyshev-Lobatto nodes of [0,1],
% ascending, and the differentiation matrix on them.
xx = cos(pi*(0:N)'/N);
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(xx,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+eye(N+1));
D = D-diag(sum(D,2));
% x = (1-xx)/2 ascends in the index order of xx, so d/dx = -2 d/dxx
x = (1-xx)/2;
D = -2*D;
end

function w = ccw(N)
% W = CCW(N) returns the Clenshaw-Curtis weights of the nodes of CHEB01.
theta = pi*(0:N)'/N;
w = zeros(N+1,1);  v = ones(N-1,1);
if mod(N,2)==0
    w(1) = 1/(N^2-1);  w(N+1) = w(1);
    for k=1:N/2-1, v = v-2*cos(2*k*theta(2:N))/(4*k^2-1); end
    v = v-cos(N*theta(2:N))/(N^2-1);
else
    w(1) = 1/N^2;  w(N+1) = w(1);
    for k=1:(N-1)/2, v = v-2*cos(2*k*theta(2:N))/(4*k^2-1); end
end
w(2:N) = 2*v/N;
w = w/2;                                    % [-1,1] -> [0,1]; symmetric
end

function pie = heatNd_pie(N,r,bc,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = HEATND_PIE(N,R,BC,DOM) PIE of the N-D heat benchmark family
%
%   u_t = sum_{i=1}^N u_{s_i s_i} + r u   on  prod_i [a_i,b_i],
%
% with one decoupled boundary condition per direction (Jagt & Peet,
% arXiv:2508.14840v4, Ex. 2 / Ex. 21 / Sec. 7.2.1 is N = 2, BC {'DD','DN'}):
%   'DD'  u = 0 at s_i = a_i and at s_i = b_i,
%   'DN'  u = 0 at a_i,  u_{s_i} = 0 at b_i,
%   'ND'  u_{s_i} = 0 at a_i,  u = 0 at b_i.
% NN is excluded: d^2/ds^2 is not invertible on it (constant mode).
%
% PIE (paper Thm. 19): state v = D^delta u = u_{s1s1...sNsN}, u = T v,
%   T = T_N ... T_1,   T_i = (d^2/ds_i^2)^{-1} on the lifted 1-D domain,
%   A = r T + A0,   A0 = sum_i prod_{j~=i} T_j   (d^2/ds_i^2 T_i = I).
%
% R IS ONLY A SHIFT OF THE RATE. Since A = A0 + r T, every LPI of Cor. 35
% at (r, k) is the LPI at (0, kappa) with kappa = r + k:
%   P*A + A*P + 2k P*T = P*A0 + A0*P + 2(r+k) P*T.
% So the benchmark is parameterized by kappa; kappa* = lambda1 (r = 0), and
% the r-invariant reach is the gap k* - k-hat = lambda1 - kappa-hat.
% MEASURED (review, 09/27/2026): the SDPs at (r, k) and (0, r + k) are
% identical in 1-D, 2-D and 3-D (bitwise at r = 12; <= 3e-17 relative at
% r = 3.7, TEST_HEATND_LPI (7)). HEATND_LPI builds from A0.
% T_i has the Green kernel (paper Cor. 4/20; L = b-a)
%   G(s,th) = h(s,th) + (s-th) [th<=s],
%   DD  h = -(s-a)(b-th)/L,   DN  h = a-s,   ND  h = th-b,
% i.e. lower cell (th<=s) h+s-th, upper cell (th>s) h; on [0,1] these are
% the paper's Ex. 3 kernels (DD: th(s-1), s(th-1); DN: -th, -s).
%
% Built directly as a 1x1 'copvar' on L2[s1..sN]: T_i commute and act in
% disjoint variables, so the kernel of T is the product of 1-D kernels,
%   params{g1..gN} = kron(C^1_{g1},...,C^N_{gN}),  ZL = ZR = {[0;1]}^N,
% (ZL(s) = kron_i ZL_i(s_i), first variable slowest, sopvar.m header), and
% each term of A is the same product with the identity {1,0,0} (ZL=ZR=0)
% in direction i. Only 2^N of the 3^N cells of T are nonzero.
%
% INPUT
% - N:   number of spatial variables, 1 <= N <= 9 (names s1..sN, sorted);
% - r:   reaction coefficient (scalar); it only shifts the rate (above),
%        so r = 0 is the canonical instance;
% - bc:  (optional) 1 x N cellstr from {'DD','DN','ND'}; default
%        {'DD','DN',...,'DN'} (the paper's example for N = 2);
% - dom: (optional) N x 2, row i = [a_i b_i]; default [0 1] rows.
%
% OUTPUT struct PIE
% - T, A:   1x1 'copvar' on L2[s1..sN];
% - A0:     1x1 'copvar', A at r = 0 (A = r T + A0), used by HEATND_LPI;
% - vars:   1 x N cellstr {'s1',...}; dom (N x 2); bc; r; N;
% - op1:    1 x N struct array of the 1-D factors (ZL, ZR, C = {mult,
%           lower, upper} coefficient matrices, K(s,th) = [1 s] C [1;th]);
% - exact:  lambda (1 x N, lowest eigenvalue of -d^2/ds_i^2 under bc{i}:
%           pi^2/L^2 for DD, pi^2/(4L^2) for DN/ND), lambda1 = sum(lambda)
%           (lowest eigenvalue of -Laplacian), kstar = lambda1 - r (exact
%           decay rate of the PDE, stable iff r < lambda1; paper App. C.1
%           extended by separation of variables), M = 1/prod(lambda) (gain
%           with ||u(t)|| <= M e^{-kstar t} ||D^delta u0||), eigfun (handle,
%           phi(S) for S n x N: product of sin/cos lowest modes), Tfac =
%           prod(-1./lambda) (T phi = Tfac phi), Afac = (r-lambda1)*Tfac;
% - desc:   one-paragraph text description.
%
% Cost: 3^N sparse cells of size 2^N x 2^N; no decision variables.
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): return A0 and document that r only
%   shifts k (kappa = r + k), which the review measured bitwise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if N<1 || N>9 || N~=round(N)
    error('heatNd_pie:N','N must be an integer in 1..9 (names s1..s9 sort).')
end
if nargin<3 || isempty(bc),     bc = [{'DD'}, repmat({'DN'},1,N-1)];    end
if ischar(bc),                  bc = {bc};                              end
if nargin<4 || isempty(dom),    dom = repmat([0 1],N,1);                end
bc = upper(reshape(bc,1,[]));
if numel(bc)~=N,    error('heatNd_pie:bc','BC must have N entries.'),   end
if ~isequal(size(dom),[N 2]) || any(dom(:,2)<=dom(:,1))
    error('heatNd_pie:dom','DOM must be N x 2 with a < b.')
end
vars = arrayfun(@(i) sprintf('s%d',i),1:N,'UniformOutput',false);

% % % 1-D factors. Coefficients C(i,j) of s^(i-1) th^(j-1).
lam = zeros(1,N);   phi = cell(1,N);
op1 = struct('ZL',cell(1,N),'ZR',[],'C',[]);
for i = 1:N
    a = dom(i,1);   b = dom(i,2);   L = b-a;
    switch bc{i}
        case 'DD'   % h = -(s-a)(b-th)/L = (ab - b s - a th + s th)/L
            Ch = [a*b, -a; -b, 1]/L;
            lam(i) = pi^2/L^2;      phi{i} = @(s) sin(pi*(s-a)/L);
        case 'DN'   % h = a - s
            Ch = [a, 0; -1, 0];
            lam(i) = pi^2/(4*L^2);  phi{i} = @(s) sin(pi*(s-a)/(2*L));
        case 'ND'   % h = th - b
            Ch = [-b, 1; 0, 0];
            lam(i) = pi^2/(4*L^2);  phi{i} = @(s) cos(pi*(s-a)/(2*L));
        otherwise
            error('heatNd_pie:bc','BC ''%s'' is not one of DD, DN, ND.',bc{i})
    end
    Cst = [0 -1; 1 0];                          % s - th, lower cell only
    op1(i).ZL = [0;1];  op1(i).ZR = [0;1];
    op1(i).C  = {zeros(2), Ch + Cst, Ch};
end
Id = struct('ZL',0,'ZR',0,'C',{{1,0,0}});      % identity in one direction

% % % T and A as tensor products of the 1-D factors.
T = copvar({kronop(num2cell(op1),vars,dom)});
A0 = [];                                        % A at r = 0
for i = 1:N
    fac = num2cell(op1);    fac{i} = Id;        % prod_{j~=i} T_j (x) I_i
    Ai = copvar({kronop(fac,vars,dom)});
    if isempty(A0),     A0 = Ai;    else,   A0 = A0 + Ai;   end
end
A = r*T + A0;

ex = struct();
ex.lambda = lam;    ex.lambda1 = sum(lam);  ex.kstar = ex.lambda1 - r;
ex.M = 1/prod(lam);
ex.eigfun = @(S) eigprod(phi,S);
ex.Tfac = prod(-1./lam);    ex.Afac = (r-ex.lambda1)*ex.Tfac;

desc = sprintf(['u_t = sum_{i=1}^%d u_{s_i s_i} + %g u on %s, BC %s. PIE (Jagt & Peet ' ...
    '2508.14840 Thm. 19): v = D^(2..2) u, T = T_%d...T_1 (Green kernels), ' ...
    'A = r T + sum_i prod_{j~=i} T_j. Exact: lambda1 = %.6f, stable iff r < lambda1, ' ...
    'decay rate k* = %.6f. The Cor. 35 LPI depends on kappa = r + k only (A = A0 + rT).'], ...
    N,r,mat2str(dom),strjoin(bc,'x'),N,ex.lambda1,ex.kstar);
pie = struct('T',T,'A',A,'A0',A0,'vars',{vars},'dom',dom,'bc',{bc},'r',r,'N',N, ...
             'op1',op1,'exact',ex,'desc',desc);
end


function P = kronop(ops,vars,dom)
% 1x1 'sopvar' on L2[vars] that is the tensor product of the 1-D operators
% ops{i} (fields ZL, ZR, C{1:3}), op i acting in vars{i} (sorted). Kernel
% of cell (g1..gN) is prod_i K^i_{gi}, hence params = kron of the 1-D
% coefficients with the first variable slowest, as ZL(s) = kron_i ZL_i.
N = numel(ops);
params = cell([3*ones(1,N),1]);
for g = 1:3^N
    c = cell(1,N);  [c{:}] = ind2sub([3*ones(1,N),1],g);    gam = [c{:}];  % dir 1 fastest
    M = 1;
    for i = 1:N,    M = kron(M,ops{i}.C{gam(i)});   end
    params{g} = sparse(M);
end
ZL = cellfun(@(o) o.ZL(:),ops,'UniformOutput',false);
ZR = cellfun(@(o) o.ZR(:),ops,'UniformOutput',false);
P = sopvar(params,struct('in',{vars},'out',{vars}),ZL,ZR, ...
           struct('in',dom,'out',dom),[1 1]);
end


function y = eigprod(phi,S)
% Product of the 1-D lowest modes, S n x N.
y = ones(size(S,1),1);
for i = 1:numel(phi),   y = y.*phi{i}(S(:,i));  end
end

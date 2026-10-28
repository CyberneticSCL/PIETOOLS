function out = opcheck_sop(mode,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = OPCHECK_SOP(MODE,...) evaluates FIXED PI operators ('copvar' or
% 'sopvar' objects) by their action on nt random polynomial test functions
% x_t: each block is applied by Gauss-Legendre quadrature of its kernel
% definition (HEATND_APPLY; sopvar document Sec. 4), and norms and inner
% products are L2 on Gauss-Legendre grids of the output spaces:
%
%   r = OPCHECK_SOP('res',L,R,OPTS)
%       max_t ||(L-R) x_t|| / max_t max(||L x_t||,||R x_t||), the relative
%       residual of the relation L = R;
%   m = OPCHECK_SOP('psd',P,OPTS)
%       min_t <x_t,P x_t> / (||x_t|| ||P x_t||), so m >= 0 exactly when
%       <x_t,P x_t> >= 0 for every t;
%   r = OPCHECK_SOP('weak',TERMS,X,OPTS)
%       the relative residual of sum_k c_k U_k' V_k = 0 in weak form, with
%       no adjoint or composition computed,
%         max_{t,u} |sum_k c_k <U_k y_u, V_k x_t>| / max_{t,u} sum_k |c_k|
%                   ||U_k y_u|| ||V_k x_t||,
%       x_t and y_u two sets of test functions on the input spaces of X;
%   c = OPCHECK_SOP('coef',P)
%       the Frobenius norm of all stored parameters of P.
%
% INPUT
% - mode:   char, 'res', 'psd', 'weak' or 'coef';
% - L, R:   'copvar' or 'sopvar' objects with the same input spaces and the
%           same output spaces;
% - P:      'copvar' or 'sopvar' object; for 'psd' its output spaces must
%           equal its input spaces;
% - TERMS:  K x 3 cell {c_1,U_1,V_1; ...; c_K,U_K,V_K}: c_k double; U_k and
%           V_k 'copvar' or 'sopvar' objects on the input spaces of X with
%           the same output spaces, or [] for the identity on the input
%           spaces of X;
% - X:      'copvar' or 'sopvar' object, whose input spaces carry x_t and
%           y_u;
% - OPTS:   (optional) 'struct' with the optional fields nt (double, the
%           number of test functions, default 6), deg (double, the largest
%           degree of a test monomial in each variable, default 3), nqi
%           (double, Gauss nodes per direction for the kernel integrals,
%           default 12), nqo (double, Gauss nodes per direction for the
%           norms and inner products, default 14) and seed (double, the rng
%           seed, default 1).
%
% OUTPUT
% - out:    double, the value r, m or c above.
%
% NOTES
% The quadrature is exact for polynomial integrands of degree at most
% 2*nqi-1 (kernel integrals) and 2*nqo-1 (norms and inner products) in each
% direction. The state of the random number generator is restored on
% return. An unknown MODE ('opcheck_sop:mode') and an operator of another
% class ('opcheck_sop:class') are errors; in 'weak', a 'sopvar' U_k with a
% nonempty V_k is not converted to a 'copvar' and is an error.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opcheck_sop
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
% Initial coding MMP, 09/29/2026. Shared by the Tier 1 end-to-end examples
%                (hinf_gain_1d_sop, stability_1d_sop, stability_2d_sop,
%                heat3d_stability_sop); generalizes the local helpers of
%                tests/test_getsol_sop (part ii) to R^n x L2 containers,
%                weak forms and a seeded generator.
% MMP, 10/21/2026: moved from examples/private/ to tests/private/ with the
%                four end-to-end examples, now the checks e2e_hinf_gain_1d,
%                e2e_stability_1d, e2e_stability_2d, e2e_heat3d_stability.
%                Code unchanged.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was a notes reference,
%                the rationale for testing against the kernel definition, a
%                cost estimate with advice on nqi and nqo, and the help
%                lines marked 10/21/2026 and their (was) line.

switch mode
    case 'res'
        % Relative residual of L = R: both on the same test functions (on
        % the input spaces of L) and the same grid (on its output spaces);
        % the generator seeded, its state restored at the end
        [L,R] = deal(as_c(varargin{1}),as_c(varargin{2}));
        o = opt(varargin,3);    s = rng;    rng(o.seed);
        F = test_funcs(L,o);    G = out_grid(L,o.nqo);
        YL = apply_c(L,F,G,o.nqi);  YR = apply_c(R,F,G,o.nqi);
        % Largest norm of (L-R) x_t over the largest norm of L x_t or R x_t
        num = 0;    den = 0;
        for t = 1:numel(F)
            num = max(num,l2n(cellfun(@minus,YL{t},YR{t},'uni',0),G));
            den = max([den,l2n(YL{t},G),l2n(YR{t},G)]);
        end
        out = num/max(den,realmin);
        rng(s);
    case 'psd'
        % P x_t on the output grid of P; the generator seeded, its state
        % restored at the end
        P = as_c(varargin{1});
        o = opt(varargin,2);    s = rng;    rng(o.seed);
        F = test_funcs(P,o);    G = out_grid(P,o.nqo);
        Y = apply_c(P,F,G,o.nqi);
        % Smallest <x_t,P x_t>/(||x_t|| ||P x_t||), x_t on the same grid,
        % as the output spaces of P equal its input spaces
        out = inf;
        for t = 1:numel(F)
            X = eval_on(F{t},G);
            out = min(out,ip(X,Y{t},G)/max(l2n(X,G)*l2n(Y{t},G),realmin));
        end
        rng(s);
    case 'weak'
        % Two sets of test functions, x_t and y_u, on the input spaces of
        % X and their grid; the generator seeded, its state restored at the
        % end
        terms = varargin{1};    X = as_c(varargin{2});
        o = opt(varargin,3);    s = rng;    rng(o.seed);
        Fx = test_funcs(X,o);   Fy = test_funcs(X,o);
        Gin = in_grid(X,o.nqo);
        nk = size(terms,1);     Ux = cell(nk,1);    Vx = cell(nk,1);    Gk = cell(nk,1);
        for k = 1:nk            % U_k y_u and V_k x_t, on the output grid of U_k (= of V_k)
            U = terms{k,2};     V = terms{k,3};
            if ~isempty(V),     V = as_c(V);    Gk{k} = out_grid(V,o.nqo);
            elseif ~isempty(U), U = as_c(U);    Gk{k} = out_grid(U,o.nqo);
            else,               Gk{k} = Gin;
            end
            Ux{k} = act(U,Fy,Gk{k},o.nqi);  Vx{k} = act(V,Fx,Gk{k},o.nqi);
        end
        % Largest |sum_k c_k <U_k y_u, V_k x_t>| over the largest
        % sum_k |c_k| ||U_k y_u|| ||V_k x_t||, over all pairs (t,u)
        num = 0;    den = 0;
        for t = 1:numel(Fx)
            for u = 1:numel(Fy)
                sv = 0;     sc = 0;
                for k = 1:nk
                    c = terms{k,1};
                    sv = sv + c*ip(Ux{k}{u},Vx{k}{t},Gk{k});
                    sc = sc + abs(c)*l2n(Ux{k}{u},Gk{k})*l2n(Vx{k}{t},Gk{k});
                end
                num = max(num,abs(sv));     den = max(den,sc);
            end
        end
        out = num/max(den,realmin);
        rng(s);
    case 'coef'
        % Frobenius norm of the parameters of every nonempty block
        P = as_c(varargin{1});  out = 0;
        for k = 1:numel(P.C)
            if isempty(P.C{k}),     continue,   end
            prm = P.C{k}.params;
            for g = 1:numel(prm),   out = out + sum(abs(prm{g}(:)).^2);    end
        end
        out = sqrt(full(out));
    otherwise
        error('opcheck_sop:mode','Unknown mode ''%s''.',mode)
end
end


% ------------------------------------------------------------------------
function o = opt(args,k)
% O = OPT(ARGS,K) returns the options: the defaults, with every field of
% the struct ARGS{K} copied over them when it is given and nonempty.
% INPUT  args: cell, the inputs after MODE; k: double, the position of OPTS
%              in ARGS.
% OUTPUT o:    'struct' with fields nt, deg, nqi, nqo and seed (double).
o = struct('nt',6,'deg',3,'nqi',12,'nqo',14,'seed',1);
if numel(args)>=k && ~isempty(args{k})
    fn = fieldnames(args{k});
    for i = 1:numel(fn),    o.(fn{i}) = args{k}.(fn{i});    end
end
end

function P = as_c(P)
% P = AS_C(P) returns a 'sopvar' block as the 1 x 1 container of itself and
% a 'copvar' unchanged; any other class is an error.
% INPUT  P: 'copvar' or 'sopvar' object.
% OUTPUT P: 'copvar' object.
if isa(P,'sopvar'),     P = copvar({P});    end
if ~isa(P,'copvar')
    error('opcheck_sop:class','Fixed operators only (copvar or sopvar); got %s.',class(P))
end
end

function F = test_funcs(P,o)
% F = TEST_FUNCS(P,O) returns O.nt random polynomial test functions on the
% input spaces of P: on input space j, 5 monomials of degree 0..O.deg in
% each variable, with N(0,1) coefficients for each of its P.dim_in(j)
% components.
% INPUT  P: 'copvar' object; o: 'struct' of options, fields nt and deg.
% OUTPUT F: 1 x o.nt cell; F{t} 1 x N cell of 'struct' with fields vars
%           (cellstr, the variables of input space j), deg (5 x numel(vars)
%           double, the exponents) and coef (5 x P.dim_in(j) double).
[~,N] = size(P.C);
F = cell(1,o.nt);
for t = 1:o.nt
    F{t} = cell(1,N);
    for j = 1:N
        vj = P.vars(P.space_in(j,:));
        F{t}{j} = struct('vars',{vj},'deg',randi([0 o.deg],5,numel(vj)),'coef',randn(5,P.dim_in(j)));
    end
end
end

function G = out_grid(P,nq)
% G = OUT_GRID(P,NQ) returns the tensor Gauss-Legendre grid of each output
% space of P.
% INPUT  P: 'copvar' object; nq: double, the nodes per direction.
% OUTPUT G: 1 x M cell; G{i} 'struct' with fields vars (cellstr), X
%           (nq^nd x nd double, the points) and w (nq^nd x 1 double, the
%           weights), nd the number of variables of output space i.
G = grid_of(P.vars,P.dom,P.space_out,nq);
end

function G = in_grid(P,nq)
% G = IN_GRID(P,NQ) returns the tensor Gauss-Legendre grid of each input
% space of P.
% INPUT  P: 'copvar' object; nq: double, the nodes per direction.
% OUTPUT G: 1 x N cell; G{j} 'struct' with fields vars (cellstr), X
%           (nq^nd x nd double, the points) and w (nq^nd x 1 double, the
%           weights), nd the number of variables of input space j.
G = grid_of(P.vars,P.dom,P.space_in,nq);
end

function G = grid_of(vars,dom,S,nq)
% G = GRID_OF(VARS,DOM,S,NQ) returns the tensor Gauss grid (points,
% weights) of each space whose variables row S(i,:) selects; an R^n space
% is one point of weight 1.
% INPUT  vars: 1 x nv cellstr, the variables; dom: nv x 2 double, their
%              domains; S: M x nv logical, row i the variables of space i;
%        nq:   double, the nodes per direction.
% OUTPUT G:    1 x M cell; G{i} 'struct' with fields vars (cellstr), X
%              (nq^nd x nd double, the points) and w (nq^nd x 1 double, the
%              weights), nd = nnz(S(i,:)).
% Nodes on [0,1]; per space, the tensor index of every point, then the
% nodes and weights mapped to the domain of each variable
[xg,wg] = gauss01(nq);
G = cell(1,size(S,1));
for i = 1:size(S,1)
    vi = vars(S(i,:));  di = dom(S(i,:),:);
    nd = numel(vi);     nx = nq^nd;
    Xp = zeros(nx,nd);  w = ones(nx,1);
    sub = cell(1,nd);
    if nd>0,    [sub{:}] = ind2sub([nq*ones(1,nd),1],(1:nx)');  end
    for t = 1:nd
        Xp(:,t) = di(t,1) + (di(t,2)-di(t,1))*xg(sub{t});
        w = w.*((di(t,2)-di(t,1))*wg(sub{t}));
    end
    G{i} = struct('vars',{vi},'X',Xp,'w',w);
end
end

function Y = act(U,F,G,nq)
% Y = ACT(U,F,G,NQ) returns U x_t on the grids G for every test function
% x_t of F; U = [] is the identity, x_t evaluated on G.
% INPUT  U:  'copvar' object or []; F: test functions (TEST_FUNCS); G:
%            grids (GRID_OF), on the output spaces of U or, for U = [], on
%            the spaces of F; nq: double, the nodes per direction of the
%            kernel integrals.
% OUTPUT Y:  1 x numel(F) cell; Y{t}{i} the values on G{i}.
if isempty(U)
    Y = cell(1,numel(F));
    for t = 1:numel(F),     Y{t} = eval_on(F{t},G);    end
else
    Y = apply_c(U,F,G,nq);
end
end

function Y = apply_c(P,F,G,nq)
% Y = APPLY_C(P,F,G,NQ) returns Y{t}{i} = (P x_t)_i on G{i}, the sum over j
% of block (i,j) applied to x_t,j by HEATND_APPLY, for every test function
% x_t.
% INPUT  P:  'copvar' object; F: test functions (TEST_FUNCS) on the input
%            spaces of P; G: grids (GRID_OF) on the output spaces of P;
%        nq: double, the nodes per direction of the kernel integrals.
% OUTPUT Y:  1 x numel(F) cell; Y{t} 1 x M cell, Y{t}{i} size(G{i}.X,1) x
%            P.dim_out(i) double.
[M,N] = size(P.C);
Y = cell(1,numel(F));
for t = 1:numel(F)
    Y{t} = cell(1,M);
    for i = 1:M
        Y{t}{i} = zeros(size(G{i}.X,1),P.dim_out(i));
        for j = 1:N
            b = P.C{i,j};
            if isempty(b),  continue,   end
            [~,ci] = ismember(b.vars.in,F{t}{j}.vars);     % columns by name
            [~,co] = ismember(b.vars.out,G{i}.vars);
            f = @(Z) eval_f(Z,ci,F{t}{j});
            Y{t}{i} = Y{t}{i} + heatNd_apply(b,f,G{i}.X(:,co),nq);
        end
    end
end
end

function X = eval_on(Ft,G)
% X = EVAL_ON(FT,G) returns the test function x_t on the grids of the same
% spaces.
% INPUT  Ft: 1 x M cell of 'struct', the components of x_t (TEST_FUNCS);
%        G:  1 x M cell of grids (GRID_OF) on the spaces of Ft.
% OUTPUT X:  1 x M cell; X{i} size(G{i}.X,1) x size(Ft{i}.coef,2) double.
X = cell(1,numel(G));
for i = 1:numel(G)
    [~,ci] = ismember(G{i}.vars,Ft{i}.vars);
    X{i} = eval_f(G{i}.X,ci,Ft{i});
end
end

function y = eval_f(Z,ci,Fj)
% Y = EVAL_F(Z,CI,FJ) returns the test-function component Fj at the points
% Z whose column k is variable Fj.vars{ci(k)}; a variable of Fj not in CI
% is zero.
% INPUT  Z:  n x numel(ci) double; ci: double, positions in Fj.vars;
%        Fj: 'struct' with fields vars, deg and coef (TEST_FUNCS).
% OUTPUT y:  n x size(Fj.coef,2) double, the sum over the monomials m of
%            prod(X.^Fj.deg(m,:)) times Fj.coef(m,:).
X = zeros(size(Z,1),numel(Fj.vars));    X(:,ci) = Z;
y = zeros(size(X,1),size(Fj.coef,2));
for m = 1:size(Fj.deg,1)
    y = y + prod(X.^Fj.deg(m,:),2).*Fj.coef(m,:);
end
end

function v = ip(X,Y,G)
% V = IP(X,Y,G) returns the L2 inner product of X and Y, summed over the
% spaces, by the quadrature weights of the grids G.
% INPUT  X, Y: 1 x M cells, X{i} and Y{i} values on G{i} of equal size;
%        G:    1 x M cell of grids (GRID_OF).
% OUTPUT v:    double.
v = 0;
for i = 1:numel(G),     v = v + sum(G{i}.w.*sum(X{i}.*Y{i},2));    end
end

function n = l2n(Y,G)
% N = L2N(Y,G) returns the L2 norm sqrt(IP(Y,Y,G)), a negative rounding of
% the square taken as 0.
% INPUT  Y: 1 x M cell of values on G; G: 1 x M cell of grids (GRID_OF).
% OUTPUT n: double.
n = sqrt(max(ip(Y,Y,G),0));
end

function [x,w] = gauss01(n)
% [X,W] = GAUSS01(N) returns the N Gauss-Legendre nodes and weights on
% [0,1] (Golub-Welsch).
% INPUT  n: double, the number of nodes.
% OUTPUT x: n x 1 double, the nodes in increasing order; w: n x 1 double,
%           the weights.
% Nodes: the eigenvalues of the Jacobi matrix of the Legendre recurrence;
% weights: twice the squared first entries of its eigenvectors; then the
% map from [-1,1] to [0,1]
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*V(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end

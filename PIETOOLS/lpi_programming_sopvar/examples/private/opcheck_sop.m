function out = opcheck_sop(mode,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = OPCHECK_SOP(MODE,...) semantic checks of FIXED operators ('copvar'
% or 'sopvar') for the examples of this folder, by their action on random
% polynomial test functions, evaluated by quadrature of the kernel
% definition (heatNd_apply, sopvar document Sec. 4). No class method acts
% on the operators here, so the checks test the operators against their
% definition, not the algebra against itself (CLAUDE.md sec. 4).
%
%   r = OPCHECK_SOP('res',L,R,OPTS)
%       max_t ||(L-R) x_t|| / max_t max(||L x_t||,||R x_t||): the relative
%       residual of the relation L = R. L and R must act on the same spaces.
%   m = OPCHECK_SOP('psd',P,OPTS)
%       min_t <x_t,P x_t> / (||x_t|| ||P x_t||) for a square P: >= 0 on
%       every test function is positivity on the test functions.
%   r = OPCHECK_SOP('weak',TERMS,X,OPTS)
%       the relation sum_k c_k U_k' V_k = 0 in weak form, with no adjoint
%       or composition computed:
%         max_{t,u} |sum_k c_k <U_k y_u, V_k x_t>| / max_{t,u} sum_k |c_k|
%                   ||U_k y_u|| ||V_k x_t||,
%       TERMS = {c_1,U_1,V_1; c_2,U_2,V_2; ...}, U_k = [] (V_k = []) the
%       identity; x_t, y_u test functions on the input spaces of the
%       operator X (any U_k or V_k), the same spaces for every term.
%   c = OPCHECK_SOP('coef',P)
%       Frobenius norm of all stored parameters of P.
%
% OPTS (struct, all optional): nt (6 test functions), deg (3, the largest
% degree per variable of a test monomial), nqi (12 Gauss nodes per
% direction for the kernel integrals), nqo (14 for the norms and inner
% products), seed (1: rng seed, restored afterwards).
%
% Quadrature: exact for polynomial integrands of degree <= 2*nq-1 per
% direction; the defaults cover kernels up to degree ~20 in 1-D. Cost per
% application: (grid points) * 3^n3 * nqi^(n3+n1) kernel evaluations per
% block, n3 shared and n1 input-only variables (heatNd_apply); use smaller
% nqi, nqo in 2-D and 3-D.
%
% Initial coding MMP, 09/29/2026. Shared by the Tier 1 end-to-end examples
%                (hinf_gain_1d_sop, stability_1d_sop, stability_2d_sop,
%                heat3d_stability_sop); generalizes the local helpers of
%                tests/test_getsol_sop (part ii) to R^n x L2 containers,
%                weak forms and a seeded generator.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch mode
    case 'res'
        [L,R] = deal(as_c(varargin{1}),as_c(varargin{2}));
        o = opt(varargin,3);    s = rng;    rng(o.seed);
        F = test_funcs(L,o);    G = out_grid(L,o.nqo);
        YL = apply_c(L,F,G,o.nqi);  YR = apply_c(R,F,G,o.nqi);
        num = 0;    den = 0;
        for t = 1:numel(F)
            num = max(num,l2n(cellfun(@minus,YL{t},YR{t},'uni',0),G));
            den = max([den,l2n(YL{t},G),l2n(YR{t},G)]);
        end
        out = num/max(den,realmin);
        rng(s);
    case 'psd'
        P = as_c(varargin{1});
        o = opt(varargin,2);    s = rng;    rng(o.seed);
        F = test_funcs(P,o);    G = out_grid(P,o.nqo);
        Y = apply_c(P,F,G,o.nqi);
        out = inf;
        for t = 1:numel(F)
            X = eval_on(F{t},G);
            out = min(out,ip(X,Y{t},G)/max(l2n(X,G)*l2n(Y{t},G),realmin));
        end
        rng(s);
    case 'weak'
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
o = struct('nt',6,'deg',3,'nqi',12,'nqo',14,'seed',1);
if numel(args)>=k && ~isempty(args{k})
    fn = fieldnames(args{k});
    for i = 1:numel(fn),    o.(fn{i}) = args{k}.(fn{i});    end
end
end

function P = as_c(P)
% A block is the 1 x 1 container of itself.
if isa(P,'sopvar'),     P = copvar({P});    end
if ~isa(P,'copvar')
    error('opcheck_sop:class','Fixed operators only (copvar or sopvar); got %s.',class(P))
end
end

function F = test_funcs(P,o)
% F{t}{j}: random polynomial on input space j of P: 5 monomials, degree
% 0..o.deg per variable, N(0,1) coefficients per component.
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
G = grid_of(P.vars,P.dom,P.space_out,nq);
end

function G = in_grid(P,nq)
G = grid_of(P.vars,P.dom,P.space_in,nq);
end

function G = grid_of(vars,dom,S,nq)
% Tensor Gauss grid (points, weights) on each space; R^n is one point.
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
% U x_t on G, U = [] the identity.
if isempty(U)
    Y = cell(1,numel(F));
    for t = 1:numel(F),     Y{t} = eval_on(F{t},G);    end
else
    Y = apply_c(U,F,G,nq);
end
end

function Y = apply_c(P,F,G,nq)
% Y{t}{i} = (P x_t)_i on G{i} = sum_j (block (i,j) applied to x_t,j).
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
% Test function x_t on the grids of the same spaces.
X = cell(1,numel(G));
for i = 1:numel(G)
    [~,ci] = ismember(G{i}.vars,Ft{i}.vars);
    X{i} = eval_f(G{i}.X,ci,Ft{i});
end
end

function y = eval_f(Z,ci,Fj)
% Fj at points Z whose column k is variable Fj.vars{ci(k)}.
X = zeros(size(Z,1),numel(Fj.vars));    X(:,ci) = Z;
y = zeros(size(X,1),size(Fj.coef,2));
for m = 1:size(Fj.deg,1)
    y = y + prod(X.^Fj.deg(m,:),2).*Fj.coef(m,:);
end
end

function v = ip(X,Y,G)
v = 0;
for i = 1:numel(G),     v = v + sum(G{i}.w.*sum(X{i}.*Y{i},2));    end
end

function n = l2n(Y,G)
n = sqrt(max(ip(Y,Y,G),0));
end

function [x,w] = gauss01(n)
% Gauss-Legendre nodes and weights on [0,1] (Golub-Welsch).
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*V(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end

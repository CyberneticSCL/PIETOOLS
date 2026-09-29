function y = heatNd_apply(P,f,S,nq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = HEATND_APPLY(P,F,S,NQ) applies a FIXED PI operator to a function by
% Gauss quadrature, straight from the 'sopvar' class definition (sopvar.m
% header; sopvar document Sec. 4), with no class method in the loop, so
% constructors and algebra are tested against the definition and not
% against each other (CLAUDE.md sec. 4):
%
%   (P x)(s) = sum_gam int K_gam(s,s') I_gam(s3-s3') x(s') ds',
%   K_gam(s,s') = (I_q (x) ZL(s))' C_gam (I_p (x) ZR(s')),
%
% gam over sorted S3 = vars.in & vars.out, direction 1 fastest; gam_k = 1
% delta (s3'_k = s3_k), 2 lower int_a^{s_k}, 3 upper int_{s_k}^b; S1
% (input only) integrated over its domain. ZL(s) = kron_i s_i.^ZL{i} over
% vars.out, ZR(s') likewise over vars.in, first variable slowest; rows of
% C_gam are (component outer, ZL monomial inner), columns likewise.
%
% INPUT
% - P:  'sopvar', or a 1x1 'copvar' (its block is used);
% - f:  handle, f(X) with X an n x numel(P.vars.in) array of points,
%       columns in P.vars.in order; returns n x p (p = P.dims(2));
% - S:  ns x numel(P.vars.out) evaluation points, columns in vars.out order;
% - nq: Gauss-Legendre nodes per direction (default 12); exact for
%       polynomial integrands of degree <= 2*nq-1 on each subinterval.
% OUTPUT
% - y:  ns x q values, q = P.dims(1).
%
% Cost: ns * 3^n3 * nq^(n3+n1) kernel evaluations; test points only.
%
% Initial coding MMP, 09/27/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(P,'copvar'),     P = P.C{1,1};   end
if nargin<4 || isempty(nq),     nq = 12;    end
vin = P.vars.in;    vout = P.vars.out;
S3 = sort(intersect(vin,vout));     S1 = setdiff(vin,vout);
n3 = numel(S3);
[~,i3in] = ismember(S3,vin);    [~,i3out] = ismember(S3,vout);
[~,i1in] = ismember(S1,vin);
q = P.dims(1);      p = P.dims(2);
[xg,wg] = gl01(nq);
ns = size(S,1);     y = zeros(ns,q);
for g = 1:3^n3
    C = full(P.params{g});
    if isempty(C) || ~any(C(:)),    continue,   end
    gam = ones(1,n3);
    if n3>0,    c = cell(1,n3);  [c{:}] = ind2sub([3*ones(1,n3),1],g);  gam = [c{1:n3}];  end
    NL = size(C,1)/q;   NR = size(C,2)/p;
    Cr = reshape(C,NL,q,NR,p);                              % (mon, comp) x (mon, comp)
    for is = 1:ns
        s = S(is,:);
        lo = [];    hi = [];    col = [];                   % integration directions
        for k = 1:n3
            if gam(k)==2,       lo(end+1) = P.dom.in(i3in(k),1);  hi(end+1) = s(i3out(k));        col(end+1) = i3in(k); %#ok<AGROW>
            elseif gam(k)==3,   lo(end+1) = s(i3out(k));          hi(end+1) = P.dom.in(i3in(k),2); col(end+1) = i3in(k); %#ok<AGROW>
            end
        end
        for k = 1:numel(S1)
            lo(end+1) = P.dom.in(i1in(k),1);    hi(end+1) = P.dom.in(i1in(k),2);    col(end+1) = i1in(k); %#ok<AGROW>
        end
        nd = numel(col);
        X = zeros(nq^nd,numel(vin));    w = ones(nq^nd,1);
        for k = 1:n3
            if gam(k)==1,   X(:,i3in(k)) = s(i3out(k));   end  % delta: s3' = s3
        end
        for t = 1:nd                                        % tensor Gauss grid
            rep = nq^(t-1);     blk = nq^(nd-t);
            xt = lo(t) + (hi(t)-lo(t))*xg;      wt = (hi(t)-lo(t))*wg;
            X(:,col(t)) = repmat(repelem(xt,rep),blk,1);
            w = w .* repmat(repelem(wt,rep),blk,1);
        end
        ZR = monvec(X,P.ZR);                                % nx x NR
        ZL = monvec(s,P.ZL);                                % 1 x NL
        fx = f(X);                                          % nx x p
        Lc = reshape(sum(Cr.*reshape(ZL(:),NL,1,1,1),1),q,NR,p);   % ZL' C
        for ii = 1:q
            acc = 0;
            for jj = 1:p
                acc = acc + sum(w.*(ZR*reshape(Lc(ii,:,jj),[],1)).*fx(:,jj));
            end
            y(is,ii) = y(is,ii) + acc;
        end
    end
end
end


function Z = monvec(X,Zc)
% Rows: kron over variables of X(:,i).^Zc{i}, first variable slowest.
n = size(X,1);  Z = ones(n,1);
for i = 1:numel(Zc)
    Zi = X(:,i).^reshape(Zc{i},1,[]);
    Z = reshape(reshape(Zi,n,[],1).*reshape(Z,n,1,[]),n,[]);   % Zi inner
end
end


function [x,w] = gl01(n)
% Gauss-Legendre nodes and weights on [0,1] (Golub-Welsch).
b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
[V,L] = eig(diag(b,1)+diag(b,-1));
[x,ix] = sort(diag(L));     w = 2*V(1,ix)'.^2;
x = (x+1)/2;    w = w/2;
end

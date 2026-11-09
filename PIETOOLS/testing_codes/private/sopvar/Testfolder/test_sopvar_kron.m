% Semantic tests of the methods kron and unkron of 'sopvar', and of the
% multiplication of a 'sopvar' with a numeric factor (mtimes), including
% operators with zero blocks in the shorthand [] or 0. For random operators
% P from L2^q[S3,S1] to L2^p[S2,S3] (random roles of s1,s2,s3), random
% matrices M, and random polynomials u, the output of each method
% applied to u is compared with the output computed from the definition of
% the operation, using only the action of P on polynomials (apply_sopvar,
% which implements the operator from its parameters) and polynomial
% arithmetic:
%   kron(M,P)*u:  block i is sum_j M(i,j)*P*u_j, u_j the j-th block of u;
%   kron(P,M)*u:  block a is sum_j P(a,:)*v_j *M(:,j),  v_j=u((0:q-1)*n+j);
%   M*P*u, P*Z*u: M*(P*u) and P*(Z*u);
%   unkron(X,n):  X*u with the copies u_i=u((0:q-1)*n+i) of the components
%                 equals interleaved copies of Xs*u_i.
% The test of unkron does not use the output of kron: X is the operator
% kron(Xs,speye(n)) with Xs random, and the output of unkron is compared
% with the action of X itself. Also tested: errors, scalar factors, the
% dimensions of the output, and that kron(M,P) and kron(P,M) for M the
% identity agree with the explicit block structure of the parameters.

% Run pietools_path_update first.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_sopvar_kron
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
% DJ (CC), 11/01/2026: Initial coding;

rng(7);
warning('off','sopvar:noncanonicalMultiplier');
names = {'s1','s2','s3'};
doms = [0,1; -1,2; 0.5,1.5];
nfail = 0;  tol = 1e-9;

% % % (1) kron(M,P) and kron(P,M), random operators (also with zero blocks)
ncase = 0;  maxerr = 0;  nshort = 0;
for t=1:200
    role = randi([0,3],1,3);
    vin = [names(role==3),names(role==1)];      vout = [names(role==2),names(role==3)];
    p = randi(3);   q = randi(3);
    [Pfull,P,ns] = rand_pair(names,doms,vin,vout,p,q,mod(t,2)==0);
    nshort = nshort+ns;
    m = randi(3);   n = randi(3);
    M = full(randn(m,n).*(rand(m,n)<0.7));
    if mod(t,3)==0, M = sparse(M); end
    % kron(M,P)
    u = rand_poly([n*q,1],vin,2);
    C = kron(M,P);
    yr = cell(m,1);
    for i=1:m
        yr{i} = polynomial(zeros(p,1));
        for j=1:n
            yr{i} = yr{i} + full(M(i,j))*apply_sopvar(Pfull,u((j-1)*q+(1:q)));
        end
    end
    e1 = pnorm(apply_sopvar(sopvar_full(C),u)-vertcat(yr{:}));
    ok = isequal(C.dims,[m*p,n*q]) && isequal(C.vars,P.vars) && isequal(C.dom.in,P.dom.in) && isequal(C.dom.out,P.dom.out);
    % kron(P,M)
    u = rand_poly([q*n,1],vin,2);
    C2 = kron(P,M);
    Y = cell(1,n);
    for j=1:n
        Y{j} = apply_sopvar(Pfull,u((0:q-1)'*n+j));
    end
    yr = cell(p,1);
    for a=1:p
        yr{a} = polynomial(zeros(m,1));
        for j=1:n
            yr{a} = yr{a} + full(M(:,j))*Y{j}(a);
        end
    end
    e2 = pnorm(apply_sopvar(sopvar_full(C2),u)-vertcat(yr{:}));
    ok = ok && isequal(C2.dims,[p*m,q*n]);
    ncase = ncase+1;    maxerr = max([maxerr,e1,e2]);
    if e1>tol || e2>tol || ~ok
        nfail = nfail+1;    fprintf('FAIL (1) t=%d role=%s p=%d q=%d m=%d n=%d e1=%.2e e2=%.2e ok=%d\n',t,mat2str(role),p,q,m,n,e1,e2,ok);
    end
end
fprintf('(1) kron(M,P), kron(P,M): %d cases (%d zero blocks in shorthand), max error %.2e\n',ncase,nshort,maxerr);

% % % (2) Numeric factors: M*P and P*Z, including zero blocks in shorthand
ncase = 0;  maxerr = 0;
for t=1:200
    role = randi([0,3],1,3);
    vin = [names(role==3),names(role==1)];      vout = [names(role==2),names(role==3)];
    p = randi(3);   q = randi(3);
    [Pfull,P] = rand_pair(names,doms,vin,vout,p,q,mod(t,2)==0);
    m = randi(3);       n = randi(3);
    M = randn(m,p);     Z = randn(q,n);
    if mod(t,5)==0,     M = sparse(M.*(rand(m,p)<0.6));     Z = sparse(Z.*(rand(q,n)<0.6));     end
    u = rand_poly([q,1],vin,2);     u2 = rand_poly([n,1],vin,2);
    QM = M*P;   QZ = P*Z;   c = randn;   QC = c*P;   QP = P*c;
    e1 = pnorm(apply_sopvar(sopvar_full(QM),u)-full(M)*apply_sopvar(Pfull,u));
    e2 = pnorm(apply_sopvar(sopvar_full(QZ),u2)-apply_sopvar(Pfull,full(Z)*u2));
    % scalar factors
    e3 = pnorm(apply_sopvar(sopvar_full(QC),u)-c*apply_sopvar(Pfull,u));
    e4 = pnorm(apply_sopvar(sopvar_full(QP),u)-c*apply_sopvar(Pfull,u));
    ok = isequal(QM.dims,[m,q]) && isequal(QZ.dims,[p,n]) && isequal(QC.dims,[p,q]) && isequal(QP.dims,[p,q]);
    ncase = ncase+1;    maxerr = max([maxerr,e1,e2,e3,e4]);
    if max([e1,e2,e3,e4])>tol || ~ok
        nfail = nfail+1;    fprintf('FAIL (2) t=%d role=%s e=%.2e %.2e %.2e %.2e ok=%d\n',t,mat2str(role),e1,e2,e3,e4,ok);
    end
end
fprintf('(2) M*P, P*Z, scalar factors: %d cases, max error %.2e\n',ncase,maxerr);

% % % (3) unkron: X=kron(Xs,speye(n)), Xs random, versus the action of X
ncase = 0;  maxerr = 0;
for t=1:100
    role = randi([0,3],1,3);
    vin = [names(role==3),names(role==1)];      vout = [names(role==2),names(role==3)];
    p = randi(3);   q = randi(3);   n = randi(3);
    [Xsf,Xs0] = rand_pair(names,doms,vin,vout,p,q,mod(t,2)==0);
    X = kron(Xs0,speye(n));
    Q = unkron(X,n);
    % action of X on u, with components (b-1)*n+i, versus Xs on the copies
    u = rand_poly([q*n,1],vin,2);
    yX = apply_sopvar(sopvar_full(X),u);
    e = 0;
    for i=1:n
        yi = apply_sopvar(sopvar_full(Q),u((0:q-1)'*n+i));
        e = max(e,pnorm(yX((0:p-1)'*n+i)-yi));
    end
    % and the extracted operator is the one of the factor
    e = max(e,pnorm(apply_sopvar(sopvar_full(Q),u(1:q))-apply_sopvar(Xsf,u(1:q))));
    ok = isequal(Q.dims,[p,q]);
    ncase = ncase+1;    maxerr = max(maxerr,e);
    if e>tol || ~ok
        nfail = nfail+1;    fprintf('FAIL (3) t=%d role=%s p=%d q=%d n=%d e=%.2e ok=%d\n',t,mat2str(role),p,q,n,e,ok);
    end
end
fprintf('(3) unkron: %d cases, max error %.2e\n',ncase,maxerr);

% % % (4) Parameter structure for the identity: kron(speye(n),P) is block
% % % diagonal in the components, kron(P,speye(n)) interleaves them
Pf = rand_op(names,doms,names([3,1]),names([2,3]),2,3);
n = 3;
A = kron(speye(n),Pf);      B = kron(Pf,speye(n));
nL = prod(cellfun(@numel,Pf.ZL));       nR = prod(cellfun(@numel,Pf.ZR));
ok = true;
for k=1:numel(Pf.params)
    Cp = Pf.params{k};
    for a=1:2
        for b=1:3
            blk = Cp((a-1)*nL+(1:nL),(b-1)*nR+(1:nR));
            for i=1:n
                ok = ok && isequal(A.params{k}(((i-1)*2+a-1)*nL+(1:nL),((i-1)*3+b-1)*nR+(1:nR)),blk);
                ok = ok && isequal(B.params{k}((((a-1)*n+i)-1)*nL+(1:nL),(((b-1)*n+i)-1)*nR+(1:nR)),blk);
            end
        end
    end
    ok = ok && nnz(A.params{k})==n*nnz(Cp) && nnz(B.params{k})==n*nnz(Cp);
end
fprintf('(4) block structure of the parameters for the identity: %s\n',pf(ok));  nfail = nfail+~ok;

% % % (5) Errors
Pf = rand_op(names,doms,names([3,1]),names([2,3]),2,3);
Xr = rand_op(names,doms,names([3,1]),names([2,3]),4,6);
errs = {@() kron(Pf,Pf), @() kron(Pf,polynomial([1;2])), @() kron(Pf,'a'), @() kron(Pf,ones(2,2,2)), ...
        @() unkron(Pf,2), @() unkron(Xr,2), @() unkron(Pf,0), @() unkron(Pf,1.5), @() unkron(Pf)};
ids = {'','','','','sopvar:unkron:notKron','sopvar:unkron:notKron','','',''};
ok = true;
for k=1:numel(errs)
    id = '-';   msg = '';
    try, errs{k}(); catch ME, msg = ME.message; id = ME.identifier; end
    if isempty(msg) || (~isempty(ids{k}) && ~strcmp(id,ids{k}))
        ok = false;     fprintf('FAIL (5) error case %d: id "%s", message "%s"\n',k,id,msg);
    end
end
% unkron of an operator that is not a Kronecker product
Xn = kron(rand_op(names,doms,names([3,1]),names([2,3]),2,2),speye(2));
Xn.params{1}(1,2) = Xn.params{1}(1,2)+1;
msg = '';   try, unkron(Xn,2); catch ME, msg = ME.identifier; end
ok = ok && strcmp(msg,'sopvar:unkron:notKron');
% unkron with N=1 returns the operator itself
Q1 = unkron(Pf,1);
ok = ok && isequal(Q1.dims,Pf.dims) && all(cellfun(@(a,b) isequal(a,b),Q1.params,Pf.params));
fprintf('(5) errors, unkron(X,1): %s\n',pf(ok));    nfail = nfail+~ok;

fprintf('\nTOTAL failures: %d\nDONE\n',nfail);

%% ========================================================================
function s = pf(ok), if ok, s = 'ok'; else, s = '<-- FAIL'; end, end

function e = pnorm(P)
P = polynomial(P);
if isempty(P.coeff), e = 0; else, e = norm(full(P.coeff),'fro'); end
end

function Q = sopvar_full(P)
% Same operator with every zero-block shorthand ([] or 0) expanded, as the
% reference action apply_sopvar reads the parameters directly
nL = prod(cellfun(@numel,P.ZL));    nR = prod(cellfun(@numel,P.ZR));
prm = P.params;
for k=1:numel(prm)
    if ~isequal(size(prm{k}),[P.dims(1)*nL,P.dims(2)*nR])
        prm{k} = sparse(P.dims(1)*nL,P.dims(2)*nR);
    end
end
Q = sopvar(prm,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end

function [Pfull,Pshort,ns] = rand_pair(names,doms,vin,vout,p,q,useshort)
% Random sopvar Pfull, and the same operator Pshort in which (if useshort)
% some zero parameters are in the shorthand [] or scalar 0; ns of them
P = rand_op(names,doms,vin,vout,p,q);
prm = P.params;     short = prm;    ns = 0;
if useshort
    for k=2:numel(prm)                              % parameter 1 stays
        if rand<0.5
            prm{k} = sparse(size(prm{k},1),size(prm{k},2));
            short{k} = [];
            if rand<0.5,    short{k} = 0;   end
            ns = ns+1;
        end
    end
end
Pfull = sopvar(prm,P.vars,P.ZL,P.ZR,P.dom,P.dims);
Pshort = sopvar(short,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end

function P = rand_op(names,doms,vin,vout,p,q)
% Random sopvar from L2^q[vin] to L2^p[vout], using rand_sopvar
[~,li] = ismember(vin,names);   [~,lo] = ismember(vout,names);
dom = struct('in',doms(li,:),'out',doms(lo,:));
P = rand_sopvar([p,q],struct('in',{vin},'out',{vout}),dom,2,0.6);
end

% Semantic tests of the 'subs' methods of sopvar and copvar. For random
% operators P with input variables S3,S1 and output variables S2,S3, all of
% the form of 'sopvar', the output of Q=subs(P,s_k,val) applied to a random
% polynomial u
% is compared with the output of P applied to u, evaluated at s_k=val using
% the 'polynomial' class (both through apply_sopvar, which implements the
% definition of the operator from its parameters, independently of 'subs').
% The multiplier parameters along s_k are made to vanish at val by projecting
% them. Also tested: the errors for non-vanishing multipliers and interior
% points, evaluation of several variables, a variable in S2, and containers.

% Run pietools_path_update first.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_subs_sopvar
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

rng(5);
warning('off','sopvar:noncanonicalMultiplier');
names = {'s1','s2','s3'};
doms = [0,1; -1,2; 0.5,1.5];
nfail = 0;  ncase = 0;  maxerr = 0;  cnt = zeros(1,4);   % S2, S3 (int a-s), S3 (int s-b), multiplier

% % % (1) Random operators, evaluated at one variable
for t=1:300
    role = randi([0,3],1,3);                 % 0: none, 1: S1, 2: S2, 3: S3
    cand = find(role>=2);
    if isempty(cand), continue, end
    k = cand(randi(numel(cand)));
    val = doms(k,1+randi([0,1]));
    vin = [names(role==3),names(role==1)];      vout = [names(role==2),names(role==3)];
    p = randi(2);   q = randi(2);
    P = rand_op(names,doms,vin,vout,p,q,names(k),val);
    u = rand_poly([q,1],vin,2);
    ye = subs(apply_sopvar(P,u),pvar(names{k}),val);
    Q = subs(P,names{k},val);
    ok = ~ismember(names{k},Q.vars.out) && isequal(sort(Q.vars.out),sort(setdiff(vout,names(k)))) && isequal(Q.dims,P.dims);
    err = pnorm(apply_sopvar(Q,u)-ye);
    ncase = ncase+1;  maxerr = max(maxerr,err);
    if err>1e-9 || ~ok
        nfail = nfail+1;  fprintf('FAIL t=%d role=%s k=%d val=%g err=%.2e ok=%d\n',t,mat2str(role),k,val,err,ok);
    end
end
fprintf('(1) random operators: %d cases, max error %.2e\n',ncase,maxerr);

% % % (2) Several variables, in a given order, and as a pvar and polynomial value
role = [3,3,2];
P = rand_op(names,doms,names([1,2]),names([3,1,2]),2,2,{'s1','s2'},[doms(1,1),doms(2,2)]);
u = rand_poly([2,1],names(role==3),2);
Q = subs(P,{'s1','s2'},[doms(1,1),doms(2,2)]);
ye = subs(subs(apply_sopvar(P,u),pvar('s1'),doms(1,1)),pvar('s2'),doms(2,2));
e = pnorm(apply_sopvar(Q,u)-ye);
Q2 = subs(P,[pvar('s1');pvar('s2')],polynomial([doms(1,1);doms(2,2)]));
e2 = pnorm(apply_sopvar(Q2,u)-ye);
Q3 = subs(P,pvar('s3'),pvar('s3'));       % x=x: unchanged
ok = e<1e-9 && e2<1e-9 && isequal(Q3.vars,P.vars) && isempty(setdiff(Q.vars.out,{'s3'}));
fprintf('(2) two variables, pvar input, identity evaluation: %s (err %.1e, %.1e)\n',pf(ok),e,e2);  nfail = nfail+~ok;

% % % (3) Errors
role = [3,2,1];
P = rand_op(names,doms,names([1,3]),names([2,1]),1,1,names(1),doms(1,2));   % multiplier vanishes at s1=b
P0 = P;
% multiplier slot that does not vanish: modify the parameter with slot 1 along s1
Pm = P;   idx = find(~cellfun(@(M) nnz(M)==0,Pm.params),1);
Pm.params{1}(1,1) = 1+rand;                      % slot (1,1,..): multiplier along s1
msg = '';  try, subs(Pm,'s1',doms(1,2)); catch ME, msg = ME.identifier; end
ok1 = strcmp(msg,'sopvar:subs:multiplier');
msg = '';  try, subs(P0,'s1',mean(doms(1,:))); catch ME, msg = ME.message; end
ok2 = ~isempty(msg);
msg = '';  try, subs(P0,'s3',doms(3,1)); catch ME, msg = ME.message; end
ok3 = ~isempty(msg);
fprintf('(3) errors: nonvanishing multiplier %s, interior point %s, input-only variable %s\n',pf(ok1),pf(ok2),pf(ok3));  nfail = nfail+~ok1+~ok2+~ok3;

% % % (4) Containers: row 1 on (s1,s2), row 2 on (s1), row 3 on (s2), columns on (s1,s3) and (s2)
val = doms(1,2);
rows = {names([1,2]),names(1),names(2)};    cols = {names([1,3]),names(2)};
C = cell(3,2);    rowdim = [2,1,2];   coldim = [2,1];
for i=1:3
    for j=1:2
        if i==2 && j==2, continue, end                                  % zero block
        C{i,j} = rand_op(names,doms,cols{j},rows{i},rowdim(i),coldim(j),{'s1'},val);
    end
end
P = copvar(C);
Q = subs(P,'s1',val);
uu = {rand_poly([P.dim_in(1),1],cols{1},2),rand_poly([P.dim_in(2),1],cols{2},2)};
e = 0;
for i=1:3
    ye = polynomial(zeros(P.dim_out(i),1));  yq = ye;
    for j=1:2
        if ~isempty(C{i,j}), ye = ye+apply_sopvar(C{i,j},uu{j}); end
        if ~isempty(Q.C{i,j}), yq = yq+apply_sopvar(Q.C{i,j},uu{j}); end
    end
    if ismember(i,[1,2]), ye = subs(ye,pvar('s1'),val); end
    e = max(e,pnorm(yq-ye));
end
ok = e<1e-9 && isempty(Q.C{2,2}) && isequal(size(Q.C),[3,2]) && isequal(Q.dim_out,P.dim_out) && isequal(Q.dim_in,P.dim_in) ...
     && isequal(Q.space_out(1,:)==1,ismember(Q.vars,{'s2'})) && isequal(Q.C{3,1}.vars,P.C{3,1}.vars);
fprintf('(4) container: %s (err %.1e)\n',pf(ok),e);  nfail = nfail+~ok;

fprintf('\nTOTAL failures: %d\nDONE\n',nfail);

%% ========================================================================
function s = pf(ok), if ok, s = 'ok'; else, s = '<-- FAIL'; end, end

function e = pnorm(P)
P = polynomial(P);
if isempty(P.coeff), e = 0; else, e = norm(full(P.coeff),'fro'); end
end

function P = rand_op(names,doms,vin,vout,p,q,kname,val)
% Random sopvar from L2^q[vin] to L2^p[vout] (rand_sopvar), whose
% multipliers along the variables in the cell kname vanish at the values in
% val
[~,li] = ismember(vin,names);   [~,lo] = ismember(vout,names);
dom = struct('in',doms(li,:),'out',doms(lo,:));
P = rand_sopvar([p,q],struct('in',{vin},'out',{vout}),dom,2,0.6);
for k=1:numel(kname)
    P = vanish_mult(P,kname{k},val(k));
end
end

function P = vanish_mult(P,name,val)
% Same operator, with the multiplier parameters along the variable name
% projected so that the multiplier vanishes at name=val: for each such
% parameter, the sum over the monomials in name of val^degree times the
% coefficient is zero. Does nothing if name is not in both the input and the
% output space of P.
vars_in = P.vars.in(:)';        vars_out = P.vars.out(:)';
S3 = intersect(vars_in,vars_out);
S3_idx = find(strcmp(S3,name));
if isempty(S3_idx)
    return
end
n3 = numel(S3);
nZL_arr = cellfun(@numel,P.ZL(:)');     nZR_arr = cellfun(@numel,P.ZR(:)');
kL = find(strcmp(vars_out,name));
% Row index: (component, monomials), the first variable slowest, so that
% the parameter reshapes to an array with the last variable fastest
shp = [fliplr(nZL_arr),P.dims(1),fliplr(nZR_arr),P.dims(2)];
pos = numel(nZL_arr)-kL+1;
w = reshape(val.^P.ZL{kL}(:),[ones(1,pos-1),nZL_arr(kL),ones(1,numel(shp)-pos)]);
params = P.params;
for lidx=1:numel(params)
    slot = cell(1,n3);      [slot{:}] = ind2sub([3*ones(1,n3),1],lidx);
    slot = cell2mat(slot);
    if slot(S3_idx)~=1
        continue
    end
    T = reshape(full(params{lidx}),shp);
    scale = max(abs(T(:)));
    T = T - w.*(sum(w.*T,pos)/sum(w.^2,pos));
    % Remove rounding errors, which subs would take for a nonzero multiplier
    % if the whole parameter is zero up to rounding
    T(abs(T)<1e-12*scale) = 0;
    params{lidx} = sparse(reshape(T,prod(nZL_arr)*P.dims(1),prod(nZR_arr)*P.dims(2)));
end
P = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end

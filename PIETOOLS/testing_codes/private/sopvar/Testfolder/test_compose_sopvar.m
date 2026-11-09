% Test of compose_sopvar. For random 'sopvar' operators X_1,...,X_nX that
% each act on a single variable (as a 3-PI operator, integrating the
% variable out, or introducing it), in a random order and with random
% dimensions, the operator Yop = compose_sopvar(...) applied to a random
% polynomial u is compared with the sequential application of the stages,
% X_nX(...(X_1 u)), computed with 'apply_sopvar' (which evaluates each stage
% from its parameters, and shares no code with compose_sopvar). Variables
% that no stage acts on pass through unchanged. Also tested: the dimensions
% and variables of Yop, the identity for an empty list of stages, and the
% errors for a mismatch of the dimensions and of the variables.
%
% Run pietools_path_update first.

% DJ (CC), 11/01/2026: Initial coding;

rng(11);
names = {'s1','s2','s3'};
doms = [0,1; -1,2; 0.5,1.5];
tol = 1e-9;
nfail = 0;      ncase = 0;      maxerr = 0;     ntype = zeros(1,4);

%% Random stages in 1 to 3 variables
for trial=1:150
    % Role of each variable: 0 (passes through), 1 (3-PI stage), 2 (stage
    % integrates the variable out), 3 (stage introduces the variable), 4 (not
    % involved)
    role = randi([0,4],1,3);
    stage_vars = names(role>=1 & role<=3);
    nX = numel(stage_vars);
    perm = randperm(nX);
    dims = randi(2,1,nX+1);                 % dims(t): input dimension of stage t
    Xops = struct('op',{},'var',{});
    for t=1:nX
        k = perm(t);                        % variable of stage t
        name = stage_vars{k};
        r = role(strcmp(names,name));
        vin = {name};   vout = {name};
        if r==2
            vout = {};
        elseif r==3
            vin = {};
        end
        [~,li] = ismember(vin,names);   [~,lo] = ismember(vout,names);
        dom = struct('in',doms(li,:),'out',doms(lo,:));
        P = rand_sopvar([dims(t+1),dims(t)],struct('in',{vin},'out',{vout}),dom,2,0.6);
        Xops(end+1) = struct('op',P,'var',name);                            %#ok<SAGROW>
        ntype(r) = ntype(r)+1;
    end
    vars_in = names(role==0 | role==1 | role==2);
    Yop = compose_sopvar(Xops,vars_in,names,doms,dims(1));

    % Reference: apply the stages one after the other
    u = rand_u(dims(1),vars_in);
    y_ref = u;
    for t=1:nX
        y_ref = apply_sopvar(Xops(t).op,y_ref);
    end
    y = apply_sopvar(Yop,u);
    err = pnorm(y-y_ref);
    vars_out = names(role==0 | role==1 | role==3);
    ok = isequal(Yop.dims,[dims(end),dims(1)]) && isequal(sort(Yop.vars.in(:)),sort(vars_in(:))) ...
            && isequal(sort(Yop.vars.out(:)),sort(vars_out(:)));
    ncase = ncase+1;    maxerr = max(maxerr,err);
    if err>tol || ~ok
        nfail = nfail+1;
        fprintf('FAIL trial %d, roles %s: err %.2e, spaces ok=%d\n',trial,mat2str(role),err,ok);
    end
end
fprintf('(1) random stages: %d cases, max error %.2e (stages: %d 3-PI, %d integrated, %d introduced)\n',ncase,maxerr,ntype(1),ntype(2),ntype(3));

%% No stages: identity
u = rand_u(2,names(1:2));
Yop = compose_sopvar(struct('op',{},'var',{}),names(1:2),names,doms,2);
err = pnorm(apply_sopvar(Yop,u)-u);
ok = err<tol && isequal(Yop.dims,[2,2]);
fprintf('(2) empty list of stages gives the identity: %s (err %.1e)\n',pf(ok),err);  nfail = nfail+~ok;

%% Errors
vin = {'s1'};   vout = {'s1'};
dom = struct('in',doms(1,:),'out',doms(1,:));
X1 = rand_sopvar([2,1],struct('in',{vin},'out',{vout}),dom,1,0.6);
X2 = rand_sopvar([1,2],struct('in',{{'s2'}},'out',{{'s2'}}),struct('in',doms(2,:),'out',doms(2,:)),1,0.6);
X3 = rand_sopvar([1,1],struct('in',{{'s2'}},'out',{{'s2'}}),struct('in',doms(2,:),'out',doms(2,:)),1,0.6);
% X3 has 1 input component, but the output of X1 has 2: the dimensions of the chain differ
bad_dim = @() compose_sopvar(struct('op',{X1,X3},'var',{'s1','s2'}),{'s1','s2'},names,doms,1);
% the variable of the stage (s2) is an input of the stage, but not in vars_in
bad_var = @() compose_sopvar(struct('op',{X1,X2},'var',{'s1','s2'}),{'s1'},names,doms,1);
ok = true;
for f = {bad_dim,bad_var}
    msg = '';
    try, f{1}(); catch ME, msg = ME.message; end
    if isempty(msg)
        ok = false;
    end
end
fprintf('(3) errors for mismatched dimensions and variables: %s\n',pf(ok));  nfail = nfail+~ok;

fprintf('\nTOTAL failures: %d\nDONE\n',nfail);

%% ========================================================================
function s = pf(ok)
if ok
    s = 'ok';
else
    s = '<-- FAIL';
end
end

function e = pnorm(P)
P = polynomial(P);
if isempty(P.coeff)
    e = 0;
else
    e = norm(full(P.coeff),'fro');
end
end

function u = rand_u(dim,vars_in)
% Random polynomial of degree 2 in the variables vars_in
if isempty(vars_in)
    u = polynomial(randn(dim,1));
else
    u = rand_poly([dim,1],vars_in,2);
end
end

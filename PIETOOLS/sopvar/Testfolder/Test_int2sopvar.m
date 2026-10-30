% Test of int2sopvar. The operator is applied with 'apply_sopvar', which
% evaluates the kernel from the definition of the class, and the result is
% compared with the integral of the kernel evaluated by the 'polynomial'
% class,
%   (Pop*x)(s) = int K(s,t) x(t) dt,
% with the dummy variables t named s_d, so that K and x can depend on the
% variable s and on its dummy variable at once. The limit [s,s] is the
% multiplier K(s) x(s), for which the kernel cannot depend on t, so only x
% is evaluated at t = s.
%
% DJ (CC), 10/29/2026: Initial coding;
% DJ (CC), 10/29/2026: A kernel depending on t is refused with the limits
%                      [s,s], where it used to be evaluated at t = s.
clear;
rng(1);
pvar s1 s2 s3
tol = 1e-10;
dom1 = [-1,2];
m = 2;      n = 3;
types_all = {'ab','as','sb','ss'};

%% Kernels of each class, between spaces with S1, S2 and S3 variables
cfgs = {struct('in',{{'s1'}},         'out',{{'s1'}}),...
        struct('in',{{'s2'}},         'out',{{'s1'}}),...
        struct('in',{{'s1','s2'}},    'out',{{'s1','s3'}}),...
        struct('in',{{'s1','s2','s3'}},'out',{{'s1','s2','s3'}}),...
        struct('in',{{'s2','s1'}},    'out',{{'s3','s1'}}),...
        struct('in',{{'s1'}},         'out',{cell(1,0)}),...
        struct('in',{cell(1,0)},      'out',{{'s1','s3'}})};
for ic=1:numel(cfgs)
    vars = cfgs{ic};
    N = numel(vars.in);     M = numel(vars.out);
    dom = struct('in',repmat(dom1,N,1),'out',repmat(dom1,M,1));
    is_S3 = ismember(vars.in,vars.out);
    vars_S1 = vars.in(~is_S3);
    dvars = strcat(vars.in,'_d');               % dummy variable of each input variable
    dvars_S1 = strcat(vars_S1,'_d');
    cfs_x = randn(n,1,2*N+2);
    x = build_poly(cfs_x,vars.in,2);
    x_d = build_poly(cfs_x,dvars,2);
    % Kernel depending on s and t: quadPoly with right variables vars.in
    cfs_Q = randn(m,n,2*(M+N)+2);
    Pk_Q = build_poly(cfs_Q,[vars.out,dvars],2);
    Q = quadPoly.polynomial2quadPoly(Pk_Q,vars.out,dvars);
    Q.nt = reshape(vars.in,1,[]);
    % Kernel in s and the variables only in vars.in: polynomial
    cfs_P = randn(m,n,2*(M+numel(vars_S1))+2);
    Kpoly = build_poly(cfs_P,[vars.out,vars_S1],2);
    Pk_P = build_poly(cfs_P,[vars.out,dvars_S1],2);
    % Constant kernel
    Kc = randn(m,n);
    kernels = {Q,Pk_Q; Kpoly,Pk_P; Kc,polynomial(Kc)};
    % Limits: default (full domain), the full domain, and random combinations
    types_cases = {[],repmat({'ab'},1,N)};
    for it=1:3
        types = repmat({'ab'},1,N);
        types(is_S3) = types_all(randi(4,1,nnz(is_S3)));
        types_cases{end+1} = types;
    end
    for ik=1:size(kernels,1)
        for it=1:numel(types_cases)
            types = types_cases{it};
            % Q depends on t of every variable, so it is no multiplier      % DJ (CC), 10/29/2026
            if ik==1 && any(strcmp(types,'ss'))                             % DJ (CC), 10/29/2026
                int_dom = int_dom_of(vars.in,types,dom1,'cell');            % DJ (CC), 10/29/2026
                assert(throws(@() int2sopvar(kernels{ik,1},vars,dom,int_dom)),...
                    "No error for a multiplier with a kernel in t.")        % DJ (CC), 10/29/2026
                continue                                                    % DJ (CC), 10/29/2026
            end                                                             % DJ (CC), 10/29/2026
            if isempty(types)
                Pop = int2sopvar(kernels{ik,1},vars,dom);
                types = repmat({'ab'},1,N);
            else
                Pop = int2sopvar(kernels{ik,1},vars,dom,int_dom_of(vars.in,types,dom1,'cell'));
                if it>2
                    % The same limits specified as 'polynomial'
                    Pop_poly = int2sopvar(kernels{ik,1},vars,dom,int_dom_of(vars.in,types,dom1,'polynomial'));
                    check_equal(apply_sopvar(Pop_poly,x),apply_sopvar(Pop,x),tol);
                end
            end
            assert(isequal(Pop.dims,[m,n]),"Operator should have the dimensions of the kernel.")
            Pexp = expected_int(kernels{ik,2}*x_d,vars.in,dvars,types,dom1);
            check_equal(apply_sopvar(Pop,x),Pexp,tol);
        end
    end
end

%% A quadPoly kernel may depend on a variable the operator lacks, with degree 0
vars = struct('in',{{'s1','s2'}},'out',{{'s1'}});
dom = struct('in',repmat(dom1,2,1),'out',dom1);
cfs_x = randn(n,1,6);
x = build_poly(cfs_x,vars.in,2);
x_d = build_poly(cfs_x,{'s1_d','s2_d'},2);
Pk = build_poly(randn(m,n,10),{'s1','s1_d','s2_d'},2);
Q = quadPoly.polynomial2quadPoly(Pk,{'s1','s9'},{'s1_d','s2_d'});
Q.nt = {'s1','s2'};
check_equal(apply_sopvar(int2sopvar(Q,vars,dom,{dom1(1),'s1';dom1(1),dom1(2)}),x),...
    expected_int(Pk*x_d,vars.in,{'s1_d','s2_d'},{'as','ab'},dom1),tol);
%% A multiplier acts pointwise, so its kernel cannot depend on t of that variable
assert(throws(@() int2sopvar(Q,vars,dom,{'s1','s1';dom1(1),dom1(2)})),"No error for a multiplier with a kernel in t.")
% It may depend on t of the other variables
Pk_s = build_poly(randn(m,n,10),{'s1','s2_d'},2);
Q_s = quadPoly.polynomial2quadPoly(Pk_s,{'s1','s9'},{'s1_d','s2_d'});
Q_s.nt = {'s1','s2'};
check_equal(apply_sopvar(int2sopvar(Q_s,vars,dom,{'s1','s1';dom1(1),dom1(2)}),x),...
    expected_int(Pk_s*x_d,vars.in,{'s1_d','s2_d'},{'ss','ab'},dom1),tol);

%% Variables and domains as for copvar: the same for the input and the output
vars_c = {'s1','s2'};       dom_c = repmat(dom1,2,1);
vars_s = struct('in',{vars_c},'out',{vars_c});
dom_s = struct('in',dom_c,'out',dom_c);
Kc = randn(m,n);
lims_c = {dom1(1),'s1';'s2','s2'};          % int_a^s1, and the multiplier along s2
cfs_x = randn(n,1,6);
x = build_poly(cfs_x,vars_c,2);
Pexp = expected_int(polynomial(Kc)*build_poly(cfs_x,{'s1_d','s2_d'},2),vars_c,{'s1_d','s2_d'},{'as','ss'},dom1);
for forms = {{vars_c,dom_c},{vars_c,dom_s},{vars_s,dom_c},{vars_s,dom_s}}
    check_equal(apply_sopvar(int2sopvar(Kc,forms{1}{1},forms{1}{2},lims_c),x),Pexp,tol);
end
Pop = int2sopvar(Kc,vars_c,dom_c);
assert(isequal(sort(Pop.vars.in),vars_c) && isequal(sort(Pop.vars.out),vars_c),"Input and output variables should both be vars.")

%% int2copvar: the same operator, as the block of a 1 x 1 copvar
Pop = int2copvar(Kc,vars_c,dom_c,lims_c);
info = verify(Pop);
assert(isa(Pop,'copvar') && isequal(size(Pop),[1,1]) && info.true==1,"Output should be a valid 1 x 1 copvar.")
check_equal(apply_sopvar(Pop.C{1,1},x),Pexp,tol);
Pop = int2copvar(Kc,vars_s,dom_s);      % without integration domain
check_equal(apply_sopvar(Pop.C{1,1},x),apply_sopvar(int2sopvar(Kc,vars_s,dom_s),x),tol);

%% Inputs that are refused
Kc = randn(m,n);
bad = {Kc,vars,dom,{dom1(1),dom1(2);dom1(1),'s2'},'variable limit for s2, which is not an output variable';
       Kc,vars,dom,[0,1;dom1],                    'limit that is not a limit of the domain';
       Kc,vars,dom,{'s1','s1','s1'},              'wrong size of the integration domain';
       Kc,vars,dom,{dom1(1),'s2';dom1(1),dom1(2)}, 'limit s2 in the row of s1';
       Kc,vars,dom,[dom1(2),dom1(1);dom1],        'limit in the wrong order';
       build_poly(randn(m,n,4),{'s3'},1),vars,dom,[], 'polynomial in a variable of neither space';
       quadPoly.polynomial2quadPoly(build_poly(randn(m,n,4),{'s9'},1),{'s9'},cell(1,0)),vars,dom,[], 'quadPoly depending on a variable that is not an output variable';
       {1},vars,dom,[],                            'kernel of the wrong class';
       Kc,struct('in',{{'s1'}}),dom,[],            'variables without a field out';
       Kc,vars,struct('in',dom1,'out',dom1),[],    'domains of the wrong size';
       Kc,vars_c,dom1,[],                          'numeric domain with one row for two variables';
       Kc,5,dom,[],                                'variables of the wrong class'};
for i=1:size(bad,1)
    assert(throws(@() int2sopvar(bad{i,1},bad{i,2},bad{i,3},bad{i,4})),"No error for %s.",bad{i,5})
end
disp('Test_int2sopvar passed.')



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_equal(Pix,Pexp,tol)
% CHECK_EQUAL(PIX,PEXP,TOL) checks that two polynomial vectors agree
%
% INPUTS
% - Pix:    m x 1 'polynomial', the output of the operator;
% - Pexp:   m x 1 'polynomial', the expected output;
% - tol:    scalar relative tolerance;
%
d = Pix-Pexp;
err = full(max([abs(d.C(:));0]));
scale = full(max([abs(Pexp.C(:));1]));
assert(err<=tol*scale,"Output differs from expected by %g (scale %g).",err,scale)
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Pexp = expected_int(integrand,vars_in,dvars,types,dom1)
% PEXP = EXPECTED_INT(INTEGRAND,VARS_IN,DVARS,TYPES,DOM1) integrates the
% integrand over the dummy variables, with the limits given by TYPES
%
% INPUTS
% - integrand:  m x 1 'polynomial' in the output variables and the dummy
%               variables;
% - vars_in:    1 x N 'cellstr', the input variables of the operator;
% - dvars:      1 x N 'cellstr', the dummy variable of each input variable;
% - types:      1 x N cell of 'char', 'ab' (a to b), 'as' (a to s), 'sb' (s
%               to b) or 'ss' (the multiplier, t = s);
% - dom1:       1 x 2 'double', the domain [a,b] of all variables;
%
% OUTPUTS
% - Pexp:       m x 1 'polynomial', the integrated polynomial;
%
Pexp = integrand;
for i=1:numel(vars_in)
    dv = polynomial(dvars(i));      sv = polynomial(vars_in(i));
    switch types{i}
        case 'ab',  Pexp = int(Pexp,dv,dom1(1),dom1(2));
        case 'as',  Pexp = int(Pexp,dv,dom1(1),sv);
        case 'sb',  Pexp = int(Pexp,dv,sv,dom1(2));
        case 'ss',  Pexp = subs(Pexp,dv,sv);
    end
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function int_dom = int_dom_of(vars_in,types,dom1,form)
% INT_DOM = INT_DOM_OF(VARS_IN,TYPES,DOM1,FORM) returns the domain of
% integration for the types of limits
%
% INPUTS
% - vars_in:    1 x N 'cellstr', the input variables of the operator;
% - types:      1 x N cell of 'char', as for 'expected_int';
% - dom1:       1 x 2 'double', the domain [a,b] of all variables;
% - form:       'char', 'cell' (variables by name) or 'polynomial';
%
% OUTPUTS
% - int_dom:    N x 2 'cell' or 'polynomial', the limits [L,U] of each
%               variable;
%
N = numel(vars_in);
lims = cell(N,2);
for i=1:N
    switch types{i}
        case 'ab',  lims(i,:) = {dom1(1),dom1(2)};
        case 'as',  lims(i,:) = {dom1(1),vars_in{i}};
        case 'sb',  lims(i,:) = {vars_in{i},dom1(2)};
        case 'ss',  lims(i,:) = {vars_in{i},vars_in{i}};
    end
end
if strcmp(form,'cell')
    int_dom = lims;
    return
end
rows = cell(N,1);
for i=1:N
    pj = cell(1,2);
    for j=1:2
        if ischar(lims{i,j})
            pj{j} = polynomial(lims(i,j));      % the variable, as pvar
        else
            pj{j} = polynomial(lims{i,j});
        end
    end
    rows{i} = [pj{:}];
end
int_dom = vertcat(rows{:});
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function P = build_poly(cfs,vnames,deg)
% P = BUILD_POLY(CFS,VNAMES,DEG) returns a polynomial matrix with given
% coefficients, so that the same coefficients can be used with other names
% for the variables
%
% INPUTS
% - cfs:    m x n x r array of 'double', the coefficients, with r at least
%           1+deg*numel(vnames)+1;
% - vnames: 1 x k cellstr, the variables of the polynomial;
% - deg:    scalar maximal degree in each variable;
%
% OUTPUTS
% - P:      m x n 'polynomial' in vnames, the sum of the constant
%           cfs(:,:,1), the powers of each variable up to deg, and the
%           product of all the variables, each with coefficients from cfs;
%
P = polynomial(cfs(:,:,1));
if isempty(vnames)
    return
end
pv = polynomial(reshape(vnames,[],1));
c = 1;
for k=1:numel(vnames)
    for d=1:deg
        c = c+1;
        P = P + cfs(:,:,c).*pv(k)^d;
    end
end
P = P + cfs(:,:,c+1).*prod(pv(:));
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function b = throws(f)
% B = THROWS(F) returns true if calling F throws an error
%
% INPUTS
% - f:      function handle without inputs;
%
% OUTPUTS
% - b:      logical, true if F() throws an error;
%
try
    f();
    b = false;
catch
    b = true;
end
end

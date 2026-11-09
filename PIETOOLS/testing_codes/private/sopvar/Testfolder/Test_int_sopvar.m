% Test of int for 'sopvar' and 'copvar'. The result is compared with the
% integral of the output of the operator, evaluated from the definition of
% the class by 'apply_sopvar' and integrated with the 'polynomial' class;
% the composition in 'int' itself is not used to determine the expected
% value.
%
% DJ (CC), 10/29/2026: Initial coding;
% DJ (CC), 10/29/2026: Added cases for the scalar U-L, when no variable of
%                      integration is an output variable of the operator.
% DJ (CC), 10/29/2026: Limits [s1,s1] now give the identity along s1, not
%                      the zero operator, so that case is changed.
clear;
rng(1);
pvar s1 s2 s3 s4
tol = 1e-10;
dom1 = [-1,2];

%% sopvar: variables, limit types, and numbers of spatial variables
% Each configuration has input and output variables, giving variables that
% are integrated out (S1), introduced (S2) and passed through (S3)
cfgs = {struct('in',{{'s1'}},         'out',{{'s1'}}),...
        struct('in',{{'s2'}},         'out',{{'s1'}}),...
        struct('in',{{'s1','s2'}},    'out',{{'s1','s3'}}),...
        struct('in',{{'s1','s2'}},    'out',{{'s1','s2'}}),...
        struct('in',{{'s1','s2','s3'}},'out',{{'s1','s2','s3'}})};
lim_types = {'ab','as','sb'};
for ic=1:numel(cfgs)
    vars = cfgs{ic};
    N = numel(vars.in);     M = numel(vars.out);
    dom = struct('in',repmat(dom1,N,1),'out',repmat(dom1,M,1));
    degs = struct('in',2*ones(1,N),'out',2*ones(1,M));
    Pop = rand_sopvar([2,3],vars,dom,degs,1);
    x = rand_poly(vars.in,3,2);
    Px = apply_sopvar(Pop,x);
    % One variable at a time, including a variable not appearing in Pop
    for v = [unique([vars.in,vars.out]),{'s4'}]
        v = v{1};   sv = polynomial({v});
        is_out = ismember(v,vars.out);
        if is_out
            lts = lim_types;
        else
            lts = {'num'};
        end
        for lt = lts
            [L,U] = lims_of(lt{1},sv,dom1);
            Pint = int(Pop,v,[L,U]);
            check_int(Pint,x,int_poly(Px,{v},[L,U]),tol);
        end
        if any(ismember(v,[vars.in,vars.out]))
            % Default domain, also for variables that are only inputs
            Pint = int(Pop,v);
            check_int(Pint,x,int_poly(Px,{v},dom1),tol);
        end
    end
    % All the output variables at once, with different limit types
    if M>0
        for it=1:3
            lt = lim_types(randi(3,1,M));
            L = cell(1,M);      U = cell(1,M);      dm = cell(M,2);
            for k=1:M
                sv = polynomial(vars.out(k));
                [L{k},U{k}] = lims_of(lt{k},sv,dom1);
                dm(k,:) = {L{k},U{k}};
            end
            Pint = int(Pop,vars.out,vertcat_poly(dm));
            check_int(Pint,x,int_poly(Px,vars.out,vertcat_poly(dm)),tol);
        end
    end
end

%% sopvar: forms of the inputs
vars = struct('in',{{'s1','s2'}},'out',{{'s1','s3'}});
dom = struct('in',repmat(dom1,2,1),'out',repmat(dom1,2,1));
Pop = rand_sopvar([2,2],vars,dom,struct('in',[2,2],'out',[2,2]),1);
x = rand_poly(vars.in,2,2);
Px = apply_sopvar(Pop,x);
Pref = int(Pop,{'s1','s3'},[dom1(1),s1;s3,dom1(2)]);
Pexp = int_poly(Px,{'s1','s3'},[dom1(1),s1;s3,dom1(2)]);
check_int(Pref,x,Pexp,tol);
forms_vars = {{'s1';'s3'},[s1;s3],[s1,s3],polynomial({'s1','s3'}),["s1","s3"]};
for i=1:numel(forms_vars)
    check_int(int(Pop,forms_vars{i},[dom1(1),s1;s3,dom1(2)]),x,Pexp,tol);
end
forms_dom = {{dom1(1),'s1';'s3',dom1(2)},{dom1(1),s1;s3,dom1(2)},{dom1(1),"s1";"s3",dom1(2)}};
for i=1:numel(forms_dom)
    check_int(int(Pop,{'s1','s3'},forms_dom{i}),x,Pexp,tol);
end
% A numeric array specifying only numbers, and no variables
check_int(int(Pop,{'s2','s4'},[0.5,1;-1,3]),x,Px*0.5*4,tol);
% Equal limits give the zero operator, also without a composition           % DJ (CC), 10/29/2026
check_int(int(Pop,'s2',[0.5,0.5]),x,0*Px,tol);                              % DJ (CC), 10/29/2026
% Integration along no variables
Pint = int(Pop,cell(1,0));
assert(isequal(Pint.vars,Pop.vars) && all(cellfun(@isequal,Pint.params,Pop.params)),'No variables should return the operator.')
%assert(all(cellfun(@nnz,Pint.params)==0),'Integral from s1 to s1 should be zero.')  % DJ (CC), 10/29/2026 (was)
% Equal variable limits, s to s, give the identity along that variable      % DJ (CC), 10/29/2026
check_int(int(Pop,'s1',[s1,s1]),x,Px,tol);                                  % DJ (CC), 10/29/2026
check_int(int(Pop,{'s1','s3'},[s1,s1;dom1(1),s3]),x,...
        int_poly(Px,{'s3'},[dom1(1),s3]),tol);                              % DJ (CC), 10/29/2026

%% sopvar: inputs that are refused
bad ={{'s3','s3'},[],'duplicate variables';
       {'s1'},[0,s3],'limit in other variable';
       {'s1'},[0,1],'limits not the domain';
       {'s2'},[0,s2],'variable limit, s2 not an output variable';
       {'s4'},[],'unknown domain, no dom';
       {'s4'},[0,s4],'variable limit, s4 not in operator';
       {'s1'},[0,1;0,1],'wrong size of dom';
       {'s1'},[0,Inf],'infinite limit';
       {'s1'},[0,s1+1],'limit in polynomial of the variable';
       {'s1'},[0,s1*s1],'limit in polynomial of the variable';
       {5},[],'variables not specified as names'};
for i=1:size(bad,1)
    if isempty(bad{i,2})
        f = @() int(Pop,bad{i,1});
    else
        f = @() int(Pop,bad{i,1},bad{i,2});
    end
    assert(throws(f),"No error for %s.",bad{i,3})
end

%% copvar: spaces with different variables, per block
spaces = struct('out',{{{},{'s1'},{'s1','s2'}}},'in',{{{'s1'},{'s1','s2'},{}}});
dims = struct('out',[2;3;1],'in',[1;2;2]);
Pop = rand_copvar(spaces,dims,dom1,1,1);
% Full domain in a variable, also if some output spaces lack the variable,
% and numeric limits in a variable not appearing in the operator. Limits
% [a,s1] and [s1,b] are only allowed if s1 appears in all output spaces
for v = {'s1','s2'}
    check_copvar(Pop,v{1},[],dom1,tol);
end
check_copvar(Pop,'s3',[0.3,1.7],dom1,tol);
check_copvar(Pop,{'s1','s2'},[],dom1,tol);
check_copvar(Pop,{'s1','s3'},[dom1;0.3,1.7],dom1,tol);
assert(throws(@() int(Pop,'s1',[dom1(1),s1])),'Variable limit in space without s1 should be refused.')
assert(throws(@() int(Pop,'s1',[0,1])),'Limits other than the domain should be refused for s1.')
spaces = struct('out',{{{'s1'},{'s1','s2'}}},'in',{{{'s1'},{'s2'},{}}});
dims = struct('out',[2;3],'in',[1;2;2]);
occ = [true,true,false;true,false,true];
Pop = rand_copvar(spaces,dims,dom1,1,1,occ);
check_copvar(Pop,'s1',[dom1(1),s1],dom1,tol);
check_copvar(Pop,'s1',[s1,dom1(2)],dom1,tol);
check_copvar(Pop,{'s1','s2'},[s1,dom1(2);dom1(1),dom1(2)],dom1,tol);
Pint = int(Pop,{'s1','s2'},[s1,dom1(2);dom1(1),dom1(2)]);
assert(all(cellfun(@isempty,Pint.C)==cellfun(@isempty,Pop.C),'all'),'Zero blocks should remain zero.')
assert(throws(@() int(Pop,'s3',[0,s3])),'Variable limit in s3 should be refused.')
% A variable only in input spaces, so in no output space: scalar b-a, which % DJ (CC), 10/29/2026
% scales Pop and keeps its metadata, and a variable limit is refused        % DJ (CC), 10/29/2026
spaces = struct('out',{{{'s1'}}},'in',{{{'s1','s2'},{}}});                  % DJ (CC), 10/29/2026
dims = struct('out',2,'in',[1;2]);                                          % DJ (CC), 10/29/2026
Pop = rand_copvar(spaces,dims,dom1,1,1);                                    % DJ (CC), 10/29/2026
check_copvar(Pop,'s2',[],dom1,tol);                                         % DJ (CC), 10/29/2026
check_copvar(Pop,{'s2','s3'},[dom1;0.3,1.7],dom1,tol);                      % DJ (CC), 10/29/2026
Pint = int(Pop,'s2');                                                       % DJ (CC), 10/29/2026
assert(isequal(Pint.vars,Pop.vars) && isequal(Pint.space_out,Pop.space_out),...
        'Metadata of Pop should be kept.')                                  % DJ (CC), 10/29/2026
assert(throws(@() int(Pop,'s2',[0,s2])),...
        'Variable limit in s2 should be refused.')                          % DJ (CC), 10/29/2026
disp('Test_int_sopvar passed.')



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_copvar(Pop,vars,dom,dom1,tol)
% CHECK_COPVAR(POP,VARS,DOM,DOM1,TOL) checks int(Pop,vars,dom) block by
% block: block (i,j) of the integral should be the integral of the output
% of block (i,j) of Pop, with respect to the variables of the output space
% i.
%
% INPUTS
% - Pop:    M x N 'copvar';
% - vars:   cellstr or char, the r variables of integration;
% - dom:    r x 2 array of limits, or [] for the default domain; a single
%           row is used for all variables;
% - dom1:   1 x 2 'double', the domain of all variables;
% - tol:    scalar tolerance;
%
vars = cellstr(vars);
dom = expand_dom(vars,dom);
Pint = int(Pop,vars,dom);
info = verify(Pint);
assert(info.true==1,"Container metadata inconsistent: %s",strjoin(info.flags,' | '))
if isempty(dom)
    dom = repmat(dom1,numel(vars),1);
end
for i=1:size(Pop.C,1)
    for j=1:size(Pop.C,2)
        assert(isempty(Pop.C{i,j})==isempty(Pint.C{i,j}),'Zero blocks should be preserved.')
        if isempty(Pop.C{i,j})
            continue
        end
        x = rand_poly(Pop.C{i,j}.vars.in,Pop.C{i,j}.dims(2),2);
        Px = apply_sopvar(Pop.C{i,j},x);
        check_int(Pint.C{i,j},x,int_poly(Px,vars,dom),tol);
    end
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function dom = expand_dom(vars,dom)
% DOM = EXPAND_DOM(VARS,DOM) repeats a single row of limits for all
% variables
%
% INPUTS
% - vars:   cellstr or char, N variables;
% - dom:    [] or 1 x 2 or N x 2 'double' or 'polynomial';
%
% OUTPUTS
% - dom:    [] or N x 2 array;
%
if ~isempty(dom) && size(dom,1)==1
    dom = repmat(dom,numel(cellstr(vars)),1);
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_int(Pint,x,Pexp,tol)
% CHECK_INT(PINT,X,PEXP,TOL) checks that the 'sopvar' PINT applied to X
% gives the expected polynomial PEXP.
%
% INPUTS
% - Pint:   'sopvar' object;
% - x:      n x 1 'polynomial' in the input variables of Pint;
% - Pexp:   m x 1 'polynomial', the expected output;
% - tol:    scalar relative tolerance;
%
Pix = apply_sopvar(Pint,x);
d = Pix-Pexp;
err = full(max([abs(d.C(:));0]));
scale = full(max([abs(Pexp.C(:));1]));
assert(err<=tol*scale,"Integral differs from expected by %g (scale %g).",err,scale)
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Pexp = int_poly(Px,vars,dom)
% PEXP = INT_POLY(PX,VARS,DOM) integrates the polynomial PX along VARS
%
% INPUTS
% - Px:     m x 1 'polynomial';
% - vars:   1 x N cellstr of the variables of integration;
% - dom:    N x 2 'double' or 'polynomial' limits [L,U], with polynomial
%           limits being the variable itself;
%
% OUTPUTS
% - Pexp:   m x 1 'polynomial', int_L^U Px dv for each variable, where
%           Px not depending on v gives (U-L)*Px;
%
Pexp = Px;
for i=1:numel(vars)
    Pexp = int(Pexp,polynomial(vars(i)),dom(i,1),dom(i,2));
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [L,U] = lims_of(type,sv,dom1)
% [L,U] = LIMS_OF(TYPE,SV,DOM1) returns the limits for the integration type
%
% INPUTS
% - type:   'ab' ([a,b]), 'as' ([a,s]), 'sb' ([s,b]), or 'num' (numbers);
% - sv:     1x1 'polynomial', the variable s;
% - dom1:   1 x 2 'double', the domain [a,b];
%
% OUTPUTS
% - L,U:    lower and upper limit, 'double' or the variable;
%
switch type
    case 'ab',  L = dom1(1);    U = dom1(2);
    case 'as',  L = dom1(1);    U = sv;
    case 'sb',  L = sv;         U = dom1(2);
    case 'num', L = 0.3;        U = 1.7;
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function dom = vertcat_poly(dm)
% DOM = VERTCAT_POLY(DM) converts an N x 2 cell of limits to a 'polynomial'
%
% INPUTS
% - dm:     N x 2 'cell' of 'double' or 'polynomial' scalars;
%
% OUTPUTS
% - dom:    N x 2 'polynomial';
%
rows = cell(size(dm,1),1);
for i=1:size(dm,1)
    rows{i} = [polynomial(dm{i,1}),polynomial(dm{i,2})];
end
dom = vertcat(rows{:});
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function x = rand_poly(vnames,n,deg)
% X = RAND_POLY(VNAMES,N,DEG) returns a random polynomial vector
%
% INPUTS
% - vnames: 1 x k cellstr, the variables of the polynomial;
% - n:      scalar number of components;
% - deg:    scalar maximal degree in each variable;
%
% OUTPUTS
% - x:      n x 1 'polynomial' in vnames, with random coefficients;
%
if isempty(vnames)
    x = polynomial(randn(n,1));
    return
end
pv = polynomial(reshape(vnames,[],1));
xc = cell(n,1);
for i=1:n
    xc{i} = polynomial(randn);
    for k=1:numel(vnames)
        for d=1:deg
            xc{i} = xc{i} + randn*pv(k)^d;
        end
    end
    % Cross term, so the variables do not separate
    xc{i} = xc{i} + randn*prod(pv(:));
end
x = vertcat(xc{:});
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

% TEST_TENSOPMAT_DIFF
% Tests @tensopmat/diff on 2 x 2 block operators
%       Cop = [ T,  M ]      T: tensopvar (pointwise product of 1-3 sopvar
%             [ K,  c ]      factors), M: polynomial matrix acting on R,
%                            K: intop (functional), c: constant matrix
% (and variants with an empty block and with a dpvar block), and on matrix-
% valued blocks M of type 'double', 'polynomial' and 'dpvar' (operators
% x -> M*x from R, function test_native_blocks) for 1-3 spatial variables.
% Oracle: apply_op(diff(Cop,vars,d),X) against the derivative of
% apply_op(Cop,X), taken by 'polynomial' calculus and stacked over the rows
% of d. The blocks are applied by apply_op/apply_sopvar/apply_op of intop,
% which share no code with the derivatives of the blocks. In
% test_native_blocks, the derivative of each block is assembled by
% 'polynomial' calculus, order by order, per row of d, which shares no code
% with the derivative of a 'dpvar'; a 'dpvar' block is built as
% M = A0 + a*A1 + b*A2 + c*A3, for polynomial matrices Ak and decision
% variables a,b,c, so that its derivative is
% D A0 + a*D A1 + b*D A2 + c*D A3.
% Cases of test_native_blocks: single blocks of random size (1x1, 2x3, 3x2,
% 1x4, 4x1), random variable lists (permuted, repeated, absent) and orders
% (including zero rows), a 2 x 2 grid of mixed native blocks, native blocks
% next to a tensopvar and a functional (compared through apply_op), input
% formats and errors.
%
% Run pietools_path_update first.

% DJ (CC), 10/20/2026: Initial coding;

rng(21);
tol = 1e-8;
npass = 0;      nfail = 0;      n_nontriv = 0;      n_err = 0;
vn = {'s1','s2'};
sv = polynomial(vn);
s1 = polynomial({'s1'});    s2 = polynomial({'s2'});
vars = struct('in',{vn},'out',{vn});
dom = [0,1];

%% Random block operators
for trial=1:25
    J = randi(3);                   % factors in T
    K = randi(2);                   % terms in T
    m1 = randi(2);      n1 = randi(2);      n2 = randi(2);
    nj = [n1;randi(2,J-1,1)];
    % T: tensopvar of J factors, all mapping L2[s1,s2]->L2[s1,s2]; the first factor
    % determines the input dimension of the block, so make all factors act on x1
    dims = [m1*ones(J,1),nj];
    depmat = struct('in',true(J,2),'out',true(J,2));
    T = rand_tensopvar(dims,vars,dom,depmat,K,[1,2],2,0.3);
    T = zero_mult_all(T,vn);
    % Block on R and constant row: the functional is tested separately below, as it acts
    % on a single function, not on a distributed monomial of J factors
    M = rand_poly([m1,n2],sv,2,5);
    c = randn(1,n2);
    Cop = tensopmat({T,M;[],c});
    % Inputs: block column 1 acts on the distributed monomial, a cell with one function
    % per factor of T, column 2 acts on R
    X1 = cell(1,J);
    for j=1:J, X1{j} = rand_poly([nj(j),1],sv,2,5); end
    x2 = randn(n2,1);
    yfun = @(Cop) apply_op(Cop,{X1,x2});
    Y = yfun(Cop);
    cases = {{{'s1'},1,       @(Y) diff(Y,s1)}, ...
             {{'s2'},1,       @(Y) diff(Y,s2)}, ...
             {{'s1','s2'},[1,1], @(Y) diff(diff(Y,s1),s2)}, ...
             {{'s2','s1'},[1,0;0,1;1,1], @(Y) [diff(Y,s2);diff(Y,s1);diff(diff(Y,s1),s2)]}, ...
             {{'s1'},2,       @(Y) diff(diff(Y,s1),s1)}, ...
             {{'s1'},[1;1],   @(Y) [diff(Y,s1);diff(Y,s1)]}};
    for c_i=1:numel(cases)
        vc = cases{c_i}{1};     dc = cases{c_i}{2};     ref_fun = cases{c_i}{3};
        try
            dCop = diff(Cop,vc,dc);
        catch ME
            % second order may fail when a factor has a nonvanishing boundary term
            if contains(ME.message,"cannot be expressed as a PI operator") && c_i==5
                n_err = n_err+1;
                continue
            else
                rethrow(ME)
            end
        end
        Yref = ref_fun(Y);
        dY = yfun(dCop);
        n_nontriv = n_nontriv + (norm_poly(Yref)>tol);
        [npass,nfail] = tally(norm_poly(dY-Yref)<tol && isequal(size(dY),size(Yref)),npass,nfail, ...
            sprintf('trial %d (J=%d,K=%d), case %d',trial,J,K,c_i));
    end
end

%% Functional block with a single-factor T
n1 = 2;     m1 = 2;     n2 = 2;
T = zero_mult_all(rand_tensopvar([m1,n1],vars,dom,struct('in',true(1,2),'out',true(1,2)),2,[1,2],2,0.3),vn);
M = rand_poly([m1,n2],sv,2,5);
Kint = rand_intop([1,n1],vn,dom,1,2,0);
c = randn(1,n2);
Cop = tensopmat({T,M;Kint,c});
x1 = rand_poly([n1,1],sv,2,5);      x2 = randn(n2,1);
Y = apply_op(Cop,{x1,x2});
dCop = diff(Cop,{'s1'},1);
Yref = [diff(Y(1:m1),s1);0*Y(m1+1)];
dY = apply_op(dCop,{x1,x2});
[npass,nfail] = tally(norm_poly(dY-Yref)<tol && isequal(size(dY),[m1+1,1]),npass,nfail,'functional and constant blocks');
% zero order returns the operator, stacked
dCop = diff(Cop,{'s1'},[0;0]);
dY = apply_op(dCop,{x1,x2});
[npass,nfail] = tally(norm_poly(dY-[Y;Y])<tol,npass,nfail,'zero order, stacked');
% variable that is not an output
dCop = diff(Cop,{'s9'},1);
[npass,nfail] = tally(norm_poly(apply_op(dCop,{x1,x2}))<tol,npass,nfail,'non-output variable');

%% Higher orders of a native block: second-order derivatives of M (the zero tensopvar block has no boundary terms)
M3 = rand_poly([m1,n2],sv,3,8);
Cz = tensopmat({0*T,M3});
Yz = apply_op(Cz,{x1,x2});
for dcase = {[2,0],[0,2],[2,1],[1,2]}
    dc = dcase{1};
    dY = apply_op(diff(Cz,{'s1','s2'},dc),{x1,x2});
    ref = Yz;
    for k=1:dc(1), ref = diff(ref,s1); end
    for k=1:dc(2), ref = diff(ref,s2); end
    [npass,nfail] = tally(norm_poly(dY-ref)<tol && norm_poly(ref)>tol,npass,nfail,sprintf('native block, order [%d,%d]',dc));
end

%% Empty block and dpvar block
Cop = tensopmat({T,[];Kint,c(:,1)});
dCop = diff(Cop,{'s1'},1);
[npass,nfail] = tally(isempty(dCop.ops{1,2}) && isequal(size(dCop.ops),[2,2]),npass,nfail,'empty block stays empty');
dpvar a b
Md = [a*s1^2 + b*s1*s2];
Td = zero_mult_all(rand_tensopvar([1,n1],vars,dom,struct('in',true(1,2),'out',true(1,2)),1,[1,2],2,0.3),vn);
Cop = tensopmat({Td,Md});
x2 = 1;
Y = apply_op(Cop,{x1,x2});
dCop = diff(Cop,{'s1'},1);
dY = apply_op(dCop,{x1,x2});
% Y = Td x1 + Md: compare the derivative of the dpvar part by its coefficients
Ypoly = apply_op(tensopmat({Td}),{x1});
dYpoly = apply_op(diff(tensopmat({Td}),{'s1'},1),{x1});
Mdref = 2*a*s1 + b*s2;
resid = dY - (dYpoly + Mdref);
[npass,nfail] = tally(isequal(size(dY),[1,1]) && max([0;abs(nonzeros(resid.C))])<tol && norm_poly(dYpoly-diff(Ypoly,s1))<tol,npass,nfail,'dpvar block');

%% Errors
Cm = tensopmat({add_mult(rand_sopvar([1,1],struct('in',{{'s1'}},'out',{{'s1'}}),[0,1],1,0.8))});
[npass,nfail] = tally(throws(@() diff(Cm,{'s1'},1),"cannot be expressed as a PI operator"),npass,nfail,'multiplier block error');
[npass,nfail] = tally(throws(@() diff(Cop,{'s1'},-1),"nonnegative"),npass,nfail,'negative order');

%% Matrix-valued native blocks
[npass,nfail] = test_native_blocks(npass,nfail);

fprintf('\nTest_tensopmat_diff: %d passed, %d failed (%d nonzero comparisons, %d expected second-order errors).\n',npass,nfail,n_nontriv,n_err);
if nfail>0
    error('Test_tensopmat_diff failed.');
end


% TEST_NATIVE_BLOCKS tests @tensopmat/diff for native blocks, see the header
function [npass,nfail] = test_native_blocks(npass,nfail)
    rng(31);
    tol = 1e-9;
    types = {'double','polynomial','dpvar'};
    sizes = {[1,1],[2,3],[3,2],[1,4],[4,1]};
    dpvar a b c

    %% Single native block, random variable lists and orders
    for nv=1:3
        vn = arrayfun(@(k) sprintf('s%d',k),1:nv,'UniformOutput',false);
        sv = polynomial(vn);
        vv = struct('in',{{}},'out',{vn});      dd = struct('in',zeros(0,2),'out',repmat([0,1],nv,1));
        for ti=1:3
            for si=1:numel(sizes)
                sz = sizes{si};
                for trial=1:4
                    B = make_block(types{ti},sz,sv);
                    Cop = tensopmat({B.M},vv,dd);
                    [vc,dc] = rand_case(vn);
                    p = size(dc,1);
                    dCop = diff(Cop,vc,dc);
                    ok = isequal(size(dCop.ops),[p,1]) && isequal(dCop.matdim,[p*sz(1),sz(2)]) && ...
                        vars_kept(dCop,Cop);
                    for k=1:p
                        blk = dCop.ops{k,1};
                        ok = ok && strcmp(class(blk),types{ti}) && blk_err(blk,exp_block(B,vc,dc(k,:)))<tol;
                    end
                    [npass,nfail] = tally(ok,npass,nfail,sprintf('single block nv=%d, %s %dx%d, trial %d, vars {%s}, d=%s',nv,types{ti},sz,trial,strjoin(vc,','),mat2str(dc)));
                end
            end
        end
    end

    %% 2 x 2 grid of native blocks of mixed type: block (i,j) of the kth derivative is block (i+(k-1)*2,j)
    for nv=1:3
        vn = arrayfun(@(k) sprintf('s%d',k),1:nv,'UniformOutput',false);
        sv = polynomial(vn);
        vv = struct('in',{{}},'out',{vn});      dd = struct('in',zeros(0,2),'out',repmat([0,1],nv,1));
        for trial=1:8
            m = randi(3,1,2);       n = randi(3,1,2);
            tp = {'dpvar','polynomial';'double','dpvar'};
            if mod(trial,2)==0, tp = tp([2,1],[2,1]); end
            B = cell(2,2);      ops = cell(2,2);
            for i=1:2
                for j=1:2
                    B{i,j} = make_block(tp{i,j},[m(i),n(j)],sv);
                    ops{i,j} = B{i,j}.M;
                end
            end
            Cop = tensopmat(ops,vv,dd);
            [vc,dc] = rand_case(vn);
            p = size(dc,1);
            dCop = diff(Cop,vc,dc);
            ok = isequal(size(dCop.ops),[2*p,2]) && isequal(dCop.matdim,[p*sum(m),sum(n)]);
            for k=1:p
                for i=1:2
                    for j=1:2
                        blk = dCop.ops{(k-1)*2+i,j};
                        ok = ok && strcmp(class(blk),tp{i,j}) && blk_err(blk,exp_block(B{i,j},vc,dc(k,:)))<tol;
                    end
                end
            end
            [npass,nfail] = tally(ok,npass,nfail,sprintf('native grid nv=%d, trial %d, vars {%s}, d=%s',nv,trial,strjoin(vc,','),mat2str(dc)));
        end
    end

    %% Native blocks next to a tensopvar (and a functional), through apply_op
    % Cop = [T, M1; [], M2],  M1: dpvar m1 x n2,  M2: polynomial m2 x n2
    for nv=1:3
        vn = arrayfun(@(k) sprintf('s%d',k),1:nv,'UniformOutput',false);
        sv = polynomial(vn);
        vars = struct('in',{vn},'out',{vn});
        vv = struct('in',{{}},'out',{vn});      dd = struct('in',zeros(0,2),'out',repmat([0,1],nv,1));
        for trial=1:4
            J = randi(2);       K = randi(2);
            m1 = randi(2);      m2 = randi(3);      n2 = randi(3);
            nj = randi(2,J,1);
            depmat = struct('in',true(J,nv),'out',true(J,nv));
            T = zero_mult_all(rand_tensopvar([m1*ones(J,1),nj],vars,[0,1],depmat,K,[1,2],2,0.3),vn);
            B1 = make_block('dpvar',[m1,n2],sv);        B2 = make_block('polynomial',[m2,n2],sv);
            Cop = tensopmat({T,B1.M;[],B2.M});
            X1 = cell(1,J);
            for j=1:J, X1{j} = rand_poly([nj(j),1],sv,2,5); end
            x2 = randn(n2,1);
            YT = apply_op(tensopmat({T}),{X1});
            [vc,dc] = rand_case(vn);
            % order 1 for the tensopvar (second order may be no PI operator)
            dc = min(dc,1);
            p = size(dc,1);
            dCop = diff(Cop,vc,dc);
            dY = apply_op(dCop,{X1,x2});
            Yref = cell(p,1);
            for k=1:p
                dYT = YT;
                for cc=1:numel(vc)
                    for l=1:dc(k,cc), dYT = diff(dYT,vc{cc}); end
                end
                Yref{k} = [dYT + exp_block(B1,vc,dc(k,:))*x2; exp_block(B2,vc,dc(k,:))*x2];
            end
            Yref = vertcat(Yref{:});
            [npass,nfail] = tally(isequal(size(dY),size(Yref)) && blk_err(dY,Yref)<tol,npass,nfail,sprintf('T + dpvar/polynomial blocks nv=%d, trial %d (J=%d,K=%d), vars {%s}, d=%s',nv,trial,J,K,strjoin(vc,','),mat2str(dc)));
        end
    end
    % Three block rows with a functional and a constant: [T, M1; [], M2; Kint, c3]
    vn = {'s1','s2'};       sv = polynomial(vn);
    vars = struct('in',{vn},'out',{vn});
    n1 = 2;     m1 = 2;     m2 = 3;     n2 = 3;
    T = zero_mult_all(rand_tensopvar([m1,n1],vars,[0,1],struct('in',true(1,2),'out',true(1,2)),2,[1,2],2,0.3),vn);
    Kint = rand_intop([1,n1],vn,[0,1],1,2,0);
    B1 = make_block('dpvar',[m1,n2],sv);    B2 = make_block('polynomial',[m2,n2],sv);    c3 = randn(1,n2);
    Cop = tensopmat({T,B1.M;[],B2.M;Kint,c3});
    x1 = rand_poly([n1,1],sv,2,5);      x2 = randn(n2,1);
    YT = apply_op(tensopmat({T}),{{x1}});
    vc = {'s2','s1'};       dc = [1,0;0,1;1,1;0,0];
    dCop = diff(Cop,vc,dc);
    dY = apply_op(dCop,{{x1},x2});
    Yref = cell(size(dc,1),1);
    for k=1:size(dc,1)
        dYT = YT;
        for cc=1:2
            for l=1:dc(k,cc), dYT = diff(dYT,vc{cc}); end
        end
        % the functional and constant rows map to R: their derivative is zero, unless
        % the order is zero (then the row is unchanged)
        if all(dc(k,:)==0)
            r3 = apply_op(tensopmat({Kint,c3}),{{x1},x2});
        else
            r3 = 0;
        end
        Yref{k} = [dYT + exp_block(B1,vc,dc(k,:))*x2; exp_block(B2,vc,dc(k,:))*x2; r3];
    end
    Yref = vertcat(Yref{:});
    [npass,nfail] = tally(isequal(size(dY),size(Yref)) && blk_err(dY,Yref)<tol,npass,nfail,'three block rows with functional and constant, stacked d');

    %% Input formats, merged variables, zero rows, errors
    vn = {'s1','s2'};       sv = polynomial(vn);        s1 = sv(1);     s2 = sv(2);
    vv = struct('in',{{}},'out',{vn});      dd = struct('in',zeros(0,2),'out',[0,1;0,1]);
    B = make_block('dpvar',[2,3],sv);       Cop = tensopmat({B.M},vv,dd);
    E1 = exp_block(B,{'s1'},1);             E12 = exp_block(B,{'s1','s2'},[1,1]);
    % pvar array, char and cellstr give the same
    d1 = diff(Cop,s1);          d2 = diff(Cop,'s1');        d3 = diff(Cop,{'s1'});       d4 = diff(Cop,'s1',1);
    [npass,nfail] = tally(all(cellfun(@(X) blk_err(X.ops{1},E1)<tol,{d1,d2,d3,d4})),npass,nfail,'input formats (default order 1)');
    d5 = diff(Cop,[s1,s2]);
    [npass,nfail] = tally(blk_err(d5.ops{1},E12)<tol,npass,nfail,'pvar array [s1,s2], default order [1,1]');
    d6 = diff(Cop,{'s2','s1'},[1,1]);
    [npass,nfail] = tally(blk_err(d6.ops{1},E12)<tol,npass,nfail,'unsorted names');
    % repeated names add their orders
    d7 = diff(Cop,{'s1','s1'},[1,2]);       d8 = diff(Cop,{'s1'},3);
    [npass,nfail] = tally(blk_err(d7.ops{1},exp_block(B,{'s1'},3))<tol && blk_err(d8.ops{1},exp_block(B,{'s1'},3))<tol,npass,nfail,'repeated names add orders');
    % scalar order for all variables
    d9 = diff(Cop,{'s1','s2'},2);
    [npass,nfail] = tally(blk_err(d9.ops{1},exp_block(B,{'s1','s2'},[2,2]))<tol,npass,nfail,'scalar order for all variables');
    % zero orders: unchanged block, in a stacked row too
    d10 = diff(Cop,{'s1','s2'},[0,0;1,0;0,0]);
    [npass,nfail] = tally(isequal(size(d10.ops),[3,1]) && blk_err(d10.ops{1},B.M)<tol && blk_err(d10.ops{3},B.M)<tol && blk_err(d10.ops{2},E1)<tol,npass,nfail,'zero rows in stacked d');
    d11 = diff(Cop,{'s1'},0);
    [npass,nfail] = tally(isequal(size(d11.ops),[1,1]) && blk_err(d11.ops{1},B.M)<tol,npass,nfail,'all-zero order returns operator');
    % variable that does not occur: zero block of the same size and class
    d12 = diff(Cop,{'s9'},1);
    [npass,nfail] = tally(isa(d12.ops{1},'dpvar') && isequal(size(d12.ops{1}),[2,3]) && blk_err(d12.ops{1},0*B.M)<tol,npass,nfail,'absent variable: zero dpvar block');
    % decision-variable-only dpvar block: zero derivative, size preserved
    Mc = [a,b,c;2*a,0*b+3,a+b];
    Cc = tensopmat({Mc},vv,dd);        dcc = diff(Cc,{'s1'},[1;2]);
    [npass,nfail] = tally(isequal(size(dcc.ops),[2,1]) && all(cellfun(@(X) isa(X,'dpvar') && isequal(size(X),[2,3]) && blk_err(X,0*Mc)<tol,dcc.ops)),npass,nfail,'decision-variable-only block');
    % output variables of a native-only row are those remaining in the derivative
    Mv = [s1*s2+s1, s1];     Cv = tensopmat({Mv},vv,dd);
    dv1 = diff(Cv,{'s1'},1);     dv2 = diff(Cv,{'s1'},2);     dv3 = diff(Cv,{'s2'},1);
    [npass,nfail] = tally(isequal(dv1.vars.out,{'s2'}) && isequal(dv1.dom.out,[0,1]) && ...
        isempty(dv2.vars.out) && isequal(size(dv2.dom.out),[0,2]) && isequal(dv3.vars.out,{'s1'}),npass,nfail,'output variables remaining in derivative');
    % errors
    [npass,nfail] = tally(throws(@() diff(Cop,{'s1'},-1),"nonnegative"),npass,nfail,'negative order');
    [npass,nfail] = tally(throws(@() diff(Cop,{'s1'},1.5),"nonnegative"),npass,nfail,'noninteger order');
    [npass,nfail] = tally(throws(@() diff(Cop,{'s1','s2'},[1,1,1]),"each of the spatial"),npass,nfail,'wrong number of orders');
    [npass,nfail] = tally(throws(@() diff(Cop,s1+s2),"polynomial combinations"),npass,nfail,'polynomial combination of variables');
    [npass,nfail] = tally(throws(@() diff(Cop,{'s1';'s2'}),"row array"),npass,nfail,'column of names');
    [npass,nfail] = tally(throws(@() diff(Cop,1),"cellstr"),npass,nfail,'numeric names');
    [npass,nfail] = tally(throws(@() diff(Cop,{'s1'},ones(2,2)),"each of the spatial"),npass,nfail,'order array with too many columns');
end


%% Helpers
function e = norm_poly(p)
if isa(p,'double')
    e = max([0;abs(p(:))]);
    return
end
p = cleanpoly(polynomial(p),1e-12);
if isempty(p.C) || nnz(p.C)==0
    e = 0;
else
    e = max(abs(p.C(:)));
end
end

function [npass,nfail] = tally(ok,npass,nfail,name)
if ok
    npass = npass+1;
else
    nfail = nfail+1;    fprintf('FAIL: %s\n',name);
end
end

function tf = throws(f,msg)
tf = false;
try
    f();
catch ME
    tf = contains(ME.message,msg);
end
end

function P = add_mult(P)
params = P.params;
params{1}(1,1) = 1;
P = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end

function Cop = zero_mult_all(Cop,vn)
% Zero all parameters of all factors with a multiplier in one of the variables vn
for l=1:numel(Cop.ops)
    P = Cop.ops{l};
    vS3 = intersect(P.vars.in,P.vars.out);
    params = P.params;
    for k=1:numel(params)
        g = gamma_of_cell(k,numel(vS3));
        if any(g(ismember(vS3,vn))==1)
            params{k} = sparse(size(params{k},1),size(params{k},2));
        end
    end
    Cop.ops{l} = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end
Cop = initialize(Cop);
end

function tf = vars_kept(dCop,Cop)
% The output variables of a row of native blocks are those remaining in its
% derivative: they must be among those of Cop, with the same domains, and
% there are no input variables
[tf,loc] = ismember(dCop.vars.out,Cop.vars.out);
tf = all(tf) && isempty(dCop.vars.in) && isequal(dCop.dom.out,Cop.dom.out(loc,:)) ...
    && isequal(size(dCop.dom.out,1),numel(dCop.vars.out));
end

function [vc,dc] = rand_case(vn)
% Random variable list (permuted, repeated, possibly absent) and orders p x r,
% p in 1:3, orders 0:3 with about 30% zeros
names = [vn,{'s9'}];
r = randi(3);
vc = names(randi(numel(names),1,r));
p = randi(3);
dc = randi([0,3],p,r);
dc(rand(p,r)<0.3) = 0;
end

function B = make_block(typ,sz,sv)
% Random block of type typ and size sz. For 'dpvar', B.A holds the polynomial
% matrices with M = A{1} + a*A{2} + b*A{3} + c*A{4}
dpvar a b c
B.typ = typ;
switch typ
    case 'double'
        B.M = randn(sz);
        B.M(rand(sz)<0.3) = 0;
    case 'polynomial'
        B.M = rand_poly(sz,sv,3,6);
    case 'dpvar'
        B.A = {rand_poly(sz,sv,3,5),rand_poly(sz,sv,2,4),rand_poly(sz,sv,2,4),rand_poly(sz,sv,3,3)};
        B.M = B.A{1} + a*B.A{2} + b*B.A{3} + c*B.A{4};
end
end

function E = exp_block(B,vc,drow)
% Expected derivative of block B, derivative of order drow(k) in variable
% vc{k} for each k, by polynomial calculus
dpvar a b c
switch B.typ
    case 'double'
        E = B.M;
        if any(drow>0), E = zeros(size(B.M)); end
    case 'polynomial'
        E = poly_diff(B.M,vc,drow);
    case 'dpvar'
        D = cellfun(@(A) poly_diff(A,vc,drow),B.A,'UniformOutput',false);
        E = D{1} + a*D{2} + b*D{3} + c*D{4};
end
end

function A = poly_diff(A,vc,drow)
for k=1:numel(vc)
    for l=1:drow(k)
        A = diff(A,vc{k});
    end
end
end

function e = blk_err(X,E)
% Largest coefficient of X-E (Inf if the sizes differ)
if ~isequal(size(X),size(E))
    e = Inf;
    return
end
if isa(X,'double') && isa(E,'double')
    e = max([0;abs(X(:)-E(:))]);
    return
end
R = X - E;
if isa(R,'dpvar')
    R = combine(R);
    e = max([0;abs(nonzeros(R.C))]);
else
    R = cleanpoly(polynomial(R),1e-12);
    e = max([0;abs(nonzeros(R.C))]);
end
end


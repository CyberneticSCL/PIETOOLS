% TEST_INTOP_DIFF
% Tests @intop/diff. A functional maps to R^m, so any derivative of
% positive order is the zero operator, with the input space of the
% functional; orders all zero return the functional, stacked p times.
% Oracle: apply_op of the result, against 0 and against the stacked
% outputs of apply_op(Kop,v), for plain 'polynomial' kernels and kernels
% with decision variables ('dpvar'), 1-3 spatial variables, varying p.
%
% Run pietools_path_update first.

% DJ (CC), 10/20/2026: Initial coding;

rng(5);
npass = 0;      nfail = 0;
tol = 1e-9;
for N=1:3
    vn = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
    sv = polynomial(vn);
    for q=[0,2]
        for trial=1:3
            m = randi(3);   n = randi(3);
            Kop = rand_intop([m,n],vn,[0,1],1,2,q);
            v = rand_poly([n,2],sv,2,5);
            Kv = apply_op(Kop,v);
            for p=1:3
                % Positive orders: zero operator
                dv = diff(Kop,vn(1),ones(p,1));
                dKv = apply_op(dv,v);
                [npass,nfail] = tally(isequal(size(dKv),[m*p,2]) && isvalid(dv) && isequal(size(dv),[m*p,n]) ...
                    && max_abs(dKv)<tol,npass,nfail,sprintf('zero, N=%d q=%d p=%d',N,q,p));
                % Several variables, higher orders
                dv = diff(Kop,vn,2*ones(p,N));
                [npass,nfail] = tally(max_abs(apply_op(dv,v))<tol && isequal(size(dv),[m*p,n]),npass,nfail,sprintf('zero, higher order, N=%d q=%d p=%d',N,q,p));
                % Zero orders: the functional itself, stacked
                dv = diff(Kop,vn(1),zeros(p,1));
                dKv = apply_op(dv,v);
                ref = repmat(Kv,p,1);
                [npass,nfail] = tally(isequal(size(dKv),[m*p,2]) && isvalid(dv) && max_abs(dKv-ref)<tol, ...
                    npass,nfail,sprintf('zero order, N=%d q=%d p=%d',N,q,p));
            end
        end
    end
end

% Input checks
Kop = rand_intop([1,2],{'s1','s2'},[0,1],1,2,0);
[npass,nfail] = tally(throws(@() diff(Kop,{'s1'},-1),"nonnegative"),npass,nfail,'negative order');
[npass,nfail] = tally(throws(@() diff(Kop,{'s1'},[1,1]),"each of the spatial"),npass,nfail,'order size');
[npass,nfail] = tally(throws(@() diff(Kop,1,1),"cellstr"),npass,nfail,'bad variable type');
% pvar input, repeated names, default order, scalar order
dv1 = diff(Kop,polynomial({'s1'}));      dv2 = diff(Kop,[polynomial({'s1'}),polynomial({'s1'})],[1,1]);
v = rand_poly([2,1],polynomial({'s1','s2'}),2,5);
[npass,nfail] = tally(max_abs(apply_op(dv1,v))<tol && max_abs(apply_op(dv2,v))<tol,npass,nfail,'pvar input, repeated names');
% a variable that is not a variable of the functional also gives zero
dv = diff(Kop,{'s9'},1);
[npass,nfail] = tally(max_abs(apply_op(dv,v))<tol,npass,nfail,'unrelated variable');

fprintf('\nTest_intop_diff: %d passed, %d failed.\n',npass,nfail);
if nfail>0
    error('Test_intop_diff failed.');
end


%% Helpers
function e = max_abs(p)
% Largest coefficient of a polynomial or dpvar array
if isa(p,'double')
    e = max(abs(p(:)));
    return
end
if isa(p,'dpvar')
    % coefficients of the (sparse) dpvar, including those of decision variables
    e = max([0;abs(nonzeros(p.C))]);
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

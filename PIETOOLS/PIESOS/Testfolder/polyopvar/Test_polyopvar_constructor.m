% Test_polyopvar_constructor.m -- verify the constructor
%   x = polyopvar(varname,dims,vars,dom,depmat)
% for vector-valued states, states depending on different subsets of the
% spatial variables, unsorted names/variables, and omitted trailing inputs.
%
% The expected field values below are written out by hand from the
% definition of each field (see the class documentation in
% @polyopvar/polyopvar.m), and the represented operator is checked by
% evaluating P(x)=C*Z(x) on test values (apply_poly) -- not by comparing
% against another PIETOOLS routine that would share the same indexing
% convention.

% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-10;
pvar s1 s2 s3

%% 1. Earlier syntax polyopvar(name,vars,dom): scalar state, full dependency
x = polyopvar('x',[s2,s1],[0,1]);
chk(isequal(x.varname,{'x'}),'1: varname');
chk(isequal(x.dims.in,1) && isequal(x.dims.out,1),'1: dims');
chk(isequal(x.vars.in,{'s1','s2'}) && isequal(x.vars.out,{'s1','s2'}),'1: vars');
chk(isequal(x.dom.in,[0,1;0,1]),'1: dom');
chk(isequal(logical(x.depmat.in),[true,true]) && isequal(logical(x.depmat.out),[true,true]),'1: depmat');

%% 2. Vector-valued states with partial dependency, unsorted names and vars
% y(s3,s2) in R^2, x(s1) in R^3, w(s3,s1) in R^1, declared in that order
vn = {'y','x','w'};
dims = [2;3;1];
vars = [s3,s1,s2];
dom = [0,3; 1,2; -1,1];                 % s3, s1, s2
depmat = logical([1 0 1; 0 1 0; 1 1 0]);    % columns: s3, s1, s2
P = polyopvar(vn,dims,vars,dom,depmat);
% States are sorted alphabetically: w, x, y
chk(isequal(P.varname,{'w';'x';'y'}),'2: varname sorted, stored as column');
chk(isequal(P.dims.in,[1;3;2]),'2: dims.in follows sorted states');
chk(isequal(P.dims.out,[2;3;1]),'2: dims.out follows rows of C (declared order)');
chk(isequal(P.vars.in,{'s1','s2','s3'}) && isequal(P.vars.out,{'s1','s2','s3'}),'2: vars');
chk(isequal(P.dom.in,[1,2;-1,1;0,3]) && isequal(P.dom.out,[1,2;-1,1;0,3]),'2: dom follows sorted vars');
% depmat.in: rows w,x,y; columns s1,s2,s3
chk(isequal(logical(P.depmat.in),logical([1 0 1; 1 0 0; 0 1 1])),'2: depmat.in');
% depmat.out: rows y,x,w (declared order); columns s1,s2,s3
chk(isequal(logical(P.depmat.out),logical([0 1 1; 1 0 0; 1 0 1])),'2: depmat.out');
chk(isequal(P.matdim,[6,1]),'2: matdim');
P2 = P;     P2.chkval = false;      P2 = initialize(P2);
chk(isequal(P2.dims,P.dims) && isequal(P2.vars,P.vars) && isequal(P2.dom,P.dom) && ...
    isequal(logical(P2.depmat.in),logical(P.depmat.in)) && isequal(logical(P2.depmat.out),logical(P.depmat.out)),...
    '2: initialize consistency');
% P(x) = [y;x;w] (declared order), for xvals ordered as P.varname
wv = rand_poly([1,1],{'s1','s3'},2,3);
xv = rand_poly([3,1],{'s1'},2,3);
yv = rand_poly([2,1],{'s2','s3'},2,3);
Y = apply_poly(P,{wv,xv,yv});
err = cleanpoly(Y - [yv;xv;wv],1e-12);
chk(isequal(size(Y),[6,1]) && max(max(abs(err.C)))<tol,'2: P(x) = [y;x;w]');

%% 3. Omitted dims (2nd input not numeric), dom shared by all variables
Q = polyopvar({'b','a'},{'s2','s1'},[0,2],logical([1 0; 0 1]));   % b(s2), a(s1)
chk(isequal(Q.varname,{'a';'b'}) && isequal(Q.dims.in,[1;1]),'3: names/dims');
chk(isequal(Q.vars.in,{'s1','s2'}) && isequal(Q.dom.in,[0,2;0,2]),'3: vars/dom');
chk(isequal(logical(Q.depmat.in),logical([1 0; 0 1])),'3: depmat.in sorted (a(s1), b(s2))');
chk(isequal(logical(Q.depmat.out),logical([0 1; 1 0])),'3: depmat.out declared order (b(s2), a(s1))');
av = rand_poly([1,1],{'s1'},2,3);      bv = rand_poly([1,1],{'s2'},2,3);
Y = apply_poly(Q,{av,bv});
err = cleanpoly(Y - [bv;av],1e-12);
chk(max(max(abs(err.C)))<tol,'3: P(x) = [b;a]');

%% 4. Trailing inputs omitted
R = polyopvar('x',3);                   % vector state, no spatial variables
chk(isequal(R.dims.in,3) && isempty(R.vars.in) && isequal(R.matdim,[3,1]),'4: dims only');
R = polyopvar({'a','b'},[2;3]);
chk(isequal(R.dims.in,[2;3]) && isequal(R.matdim,[5,1]),'4: dims only, 2 states');
R = polyopvar('x',2,{'s1'});            % dom default [0,1]
chk(isequal(R.dom.in,[0,1]) && isequal(R.dims.in,2),'4: dims, vars');
R = polyopvar('x',2,[s2,s1],[0,1;1,2]); % dom given per variable (s2,s1), sorted to (s1,s2)
chk(isequal(R.dom.in,[1,2;0,1]) && isequal(R.vars.in,{'s1','s2'}),'4: dims, vars, dom');

%% 5. vartab must keep sizes and spatial dependency of the states
V = x.vartab;
chk(isequal(V.vars.in,x.vars.in) && isequal(V.dom.in,x.dom.in),'5: vartab');
V = P.vartab;       % states of different sizes and dependency (from 2)
% vartab declares the states in sorted order (w,x,y): dims.out, depmat.out follow it, not (2)
chk(isequal(V.varname,P.varname) && isequal(V.dims.in,[1;3;2]) && isequal(V.dims.out,[1;3;2]) && ...
    isequal(V.vars,P.vars) && isequal(V.dom,P.dom),'5: vartab, mixed states: names/dims/vars/dom');
chk(isequal(logical(V.depmat.in),logical([1 0 1; 1 0 0; 0 1 1])) && ...
    isequal(logical(V.depmat.out),logical([1 0 1; 1 0 0; 0 1 1])),'5: vartab, mixed states: depmat');

%% 6. Invalid inputs
chk_error(@() polyopvar({'a','b'},[1;2;3],s1),'6: wrong number of dims');
chk_error(@() polyopvar({'a','b'},[1;1.5],s1),'6: non-integer dims');
chk_error(@() polyopvar({'a','b'},[1;0],s1),'6: zero dims');
chk_error(@() polyopvar({'a','b'},[1;1],[s1,s2],[0,1],true(2,3)),'6: wrong depmat size');
chk_error(@() polyopvar({'a','b'},[1;1],[s1,s2],[0,1;0,1;0,1]),'6: wrong dom size');

disp('Test_polyopvar_constructor passed.')

function chk(cond,msg)
    if ~cond
        error("Test_polyopvar_constructor failed: %s.",msg)
    end
end
function chk_error(fn,msg)
    try
        fn();
    catch
        return
    end
    error("Test_polyopvar_constructor failed: no error thrown (%s).",msg)
end

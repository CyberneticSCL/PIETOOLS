% Test_polyopvar_initialize_vars.m -- verify @polyopvar/initialize.m when
% f.vars.in and f.C.vars.in differ.
%
%   Supported: variable of f that C lacks (no monomial may depend on it),
%   variable of C that f lacks (C may not act on it), and reordered
%   f.vars.in. Rejected: monomial depending on a variable C lacks; C acting
%   on a variable f lacks; wrong domains.
% Oracle: expected values are set by hand from the declarations; the
% polynomial is evaluated with apply_poly before/after.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;
pvar s1 s2 s3

% Two states on different variables: x(s1), y(s1,s2), domains [0,1], [0,2]
x = polyopvar('x',[s1,s2],[0,1;0,2]);
y = polyopvar('y',[s1,s2],[0,1;0,2]);
f = x + y + x.*y;
if ~isequal(f.vars.in,{'s1','s2'}) || ~isequal(f.dom.in,[0,1;0,2]) || ...
        ~isequal(logical(f.depmat.in),true(2,2))
    error("Test_polyopvar_initialize_vars failed: unexpected reference object.")
end
xv = rand_poly([1,1],{'s1','s2'},2,5);
yv = rand_poly([1,1],{'s1','s2'},2,5);
Fref = apply_poly(f,{xv,yv});

% (1) f.vars.in reordered, dom.in and depmat.in left empty: inferred values
%     must follow the order of f.vars.in
g = f;      g.vars.in = {'s2','s1'};
g.dom.in = zeros(0,2);          g.depmat.in = false(0,2);       g.chkval = false;
g = initialize(g);
if ~isequal(g.dom.in,[0,2;0,1]) || ~isequal(logical(g.depmat.in),true(2,2))
    error("Test_polyopvar_initialize_vars failed: dom/depmat not in order of f.vars.in.")
end

% (2) f has an extra variable s3 that C lacks, no state depends on it
g = f;      g.vars.in = {'s1','s2','s3'};
g.dom.in = [0,1;0,2;0,3];       g.depmat.in = [f.depmat.in,false(2,1)];
g.chkval = false;
h = initialize(g);
if ~h.chkval || ~isequal(h.vars.in,{'s1','s2','s3'}) || ~isequal(h.dom.in,[0,1;0,2;0,3])
    error("Test_polyopvar_initialize_vars failed: extra variable of f not accepted.")
end
% ... and with depmat.in empty (inferred as no dependence on s3)
g.depmat.in = false(0,3);       g.chkval = false;
h = initialize(g);
if ~isequal(logical(h.depmat.in),logical([1,1,0;1,1,0]))
    error("Test_polyopvar_initialize_vars failed: depmat.in not inferred for extra variable.")
end
Fh = apply_poly(h,{xv,yv});
err = cleanpoly(Fh-Fref,1e-10);
if max(max(abs(err.C)))>tol*max(1,max(max(abs(Fref.C))))
    error("Test_polyopvar_initialize_vars failed: polynomial changed.")
end

% (3) Extra variable of f, but domain not given: cannot be inferred
g.depmat.in = [f.depmat.in,false(2,1)];     g.dom.in = zeros(0,2);
g.chkval = false;
expect_error(@() initialize(g),'must be specified','extra variable without domain');

% (4) Monomial depends on a variable that C lacks
g = f;      g.vars.in = {'s1','s2','s3'};
g.dom.in = [0,1;0,2;0,3];       g.depmat.in = [f.depmat.in,[1;0]];
g.chkval = false;
expect_error(@() initialize(g),'not defined','monomial depending on variable missing in C');

% (5) C has an extra variable s3 that f lacks
C = f.C;
C.vars.in = {'s1','s2','s3'};   C.dom.in = [0,1;0,2;0,3];
% 5a: C does not act on s3
C.depmat.in = [f.C.depmat.in,zeros(size(f.C.depmat.in,1),1)];
g = f;      g.C = C;            g.chkval = false;
h = initialize(g);              % f.vars.in = {s1,s2}
if ~isequal(h.vars.in,{'s1','s2'}) || ~isequal(h.dom.in,[0,1;0,2]) || ...
        ~isequal(logical(h.depmat.in),true(2,2))
    error("Test_polyopvar_initialize_vars failed: extra variable of C not accepted.")
end
g.vars.in = cell(1,0);  g.dom.in = zeros(0,2);  g.depmat.in = false(0,2);
g.chkval = false;
h = initialize(g);              % f.vars.in empty: copied from C
if ~isequal(h.vars.in,{'s1','s2','s3'})
    error("Test_polyopvar_initialize_vars failed: empty vars.in not copied from C.")
end
% 5b: C acts on s3
C.depmat.in = [f.C.depmat.in,ones(size(f.C.depmat.in,1),1)];
g = f;      g.C = C;            g.chkval = false;
expect_error(@() initialize(g),'conflicts','C acting on variable missing in f');

% (6) Wrong domain of a shared variable
g = f;      g.dom.in = [0,1;0,5];       g.chkval = false;
expect_error(@() initialize(g),'domains','wrong domain');

disp("Test_polyopvar_initialize_vars passed.")


function expect_error(fn,msgpart,label)
try
    fn();
catch ME
    if ~contains(ME.message,msgpart)
        error("Test_polyopvar_initialize_vars failed (%s): rejected with unexpected message: %s",label,ME.message)
    end
    return
end
error("Test_polyopvar_initialize_vars failed (%s): not rejected.",label)
end


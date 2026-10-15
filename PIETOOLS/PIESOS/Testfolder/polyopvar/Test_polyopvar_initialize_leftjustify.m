% Test_polyopvar_initialize_leftjustify.m -- verify @polyopvar/initialize.m
% moves the nonzeros of each row of Z to the leading columns (order kept).
%
% Oracle: Z is scrambled by hand and checked structurally, and
% independently by evaluating the polynomial (apply_poly) before/after.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
y = polyopvar('y',[s1,s2],[0,1]);
f = x + x.*y + x.*y.*x + y.*y;
Zref = f.Z;
xv = {rand_poly([1,1],{'s1','s2'},2,5)};
yv = {rand_poly([1,1],{'s1','s2'},2,5)};
Fref = apply_poly(f,{xv{1},yv{1}});

% Scramble: shift each row's nonzeros to arbitrary columns, keeping order
Zfull = full(Zref);      d = size(Zfull,2);
Zs = zeros(size(Zfull,1),d+2);
for i=1:size(Zfull,1)
    nz = Zfull(i,Zfull(i,:)~=0);
    cols = sort(randperm(d+2,numel(nz)));
    Zs(i,cols) = nz;
end
g = f;      g.Z = sparse(Zs);       g.chkval = false;
g = initialize(g);

Zg = full(g.Z);
if ~isequal(Zg(:,1:d),Zfull) || any(any(Zg(:,d+1:end)))
    error("Test_polyopvar_initialize_leftjustify failed: Z not left-justified.")
elseif ~issparse(g.Z)
    error("Test_polyopvar_initialize_leftjustify failed: Z not sparse.")
end
Fg = apply_poly(g,{xv{1},yv{1}});
err = cleanpoly(Fg-Fref,1e-10);
if max(max(abs(err.C)))>tol*max(1,max(max(abs(Fref.C))))
    error("Test_polyopvar_initialize_leftjustify failed: polynomial changed.")
end

% Already left-justified input is unchanged
h = f;      h.chkval = false;       h = initialize(h);
if ~isequal(h.Z,f.Z)
    error("Test_polyopvar_initialize_leftjustify failed: justified Z altered.")
end

disp("Test_polyopvar_initialize_leftjustify passed.")

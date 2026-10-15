% Test_polyopvar_initialize.m -- verify @polyopvar/initialize.m
% correctly infers missing fields (dims/vars/dom/depmat) from the
% coefficient operator, across varying numbers of spatial variables and
% state dimensions, and correctly rejects an internally inconsistent
% object.
%
% Oracle: the reference object is built by the constructor (x.*x,
% x+x.*x) before its fields are stripped/corrupted; values are compared
% with apply_poly.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - Negative checks test ME.message, not ME.identifier: a plain error()
%   leaves the identifier empty, so an identifier test can never fire.
% - Case B checks a known limitation: dims.in(i)>1 cannot be recovered
%   for a state that only appears at degree >= 2 (a monomial fixes only
%   the product of its factor sizes), so initialize must error. The test
%   fails if initialize recovers it.
% - Vector-valued states in general: see Test_polyopvar_vector_states.m.

tol = 1e-8;

for Ntot = 1:3
    for trial = 1:8
        rng(100*Ntot+trial);

        % Declare Ntot spatial variables and a scalar state depending on
        % a random nonempty subset of them.
        svarname = cell(1,Ntot);
        for j=1:Ntot
            svarname{j} = ['s',num2str(j)];
        end
        isvar = false(1,Ntot);
        while ~any(isvar)
            isvar = randi([0,1],1,Ntot)==1;
        end
        svars = svarname(isvar);
        dom = [zeros(numel(svars),1),ones(numel(svars),1)];

        pvarlist = cell(1,Ntot);
        for j=1:Ntot
            pvarlist{j} = pvar(svarname{j});
        end
        pvars_full = [pvarlist{:}];

        x = polyopvar('x',svars,dom);

        % Build x2 = x.*x by hand, then strip dims/vars/dom/depmat.in and
        % confirm initialize.m re-derives them identically.
        x2 = x.*x;
        x2_stripped = x2;
        x2_stripped.dims.in = zeros(1,0);
        x2_stripped.vars.in = cell(1,0);
        x2_stripped.dom.in = zeros(0,2);
        x2_stripped.depmat.in = false(1,0);
        x2_stripped.chkval = false;

        x2_init = initialize(x2_stripped);
        if ~x2_init.chkval
            error("Test_polyopvar_initialize failed: chkval not set (Ntot=%d, trial=%d).",Ntot,trial)
        end
        if ~isequal(x2_init.dims.in,x2.dims.in)
            error("Test_polyopvar_initialize failed: dims.in not correctly inferred (Ntot=%d, trial=%d).",Ntot,trial)
        end
        if ~isequal(logical(x2_init.depmat.in),logical(x2.depmat.in))
            error("Test_polyopvar_initialize failed: depmat.in not correctly inferred (Ntot=%d, trial=%d).",Ntot,trial)
        end

        % Numeric check: initialize.m must not have altered the polynomial
        xv = rand_poly([1,1],svars,2,5);
        Y_orig = apply_poly(x2,{xv});
        Y_init = apply_poly(x2_init,{xv});
        err = cleanpoly(Y_orig-Y_init,1e-10);
        if max(max(abs(err.C)))>tol
            error("Test_polyopvar_initialize failed: value changed after initialize (Ntot=%d, trial=%d).",Ntot,trial)
        end

        % Negative check: a deliberately-inconsistent dims.in must be rejected
        x2_bad = x2;
        x2_bad.dims.in = x2.dims.in+1;
        x2_bad.chkval = false;
        try
            initialize(x2_bad);
            error("Test_polyopvar_initialize failed: inconsistent dims.in was not rejected (Ntot=%d, trial=%d).",Ntot,trial)
        catch ME
            if contains(ME.message,'Test_polyopvar_initialize')
                rethrow(ME);
            end
        end
    end
end

disp('Test_polyopvar_initialize passed.')

% Vector-valued states (dims.in(i)>1): rand_polyopvar builds one state of
% random size in {1,2,3}.
for trial = 1:20
    rng(9100+trial);

    [x,~,xdims,xvars] = rand_polyopvar(1,3);
    d = xdims(1);

    % Case A: state appears bare (degree 1) in Z, so dims.in is inferred from
    % the degree-1 block of C, as in the scalar case.
    f = x + x.*x;
    f_stripped = f;
    f_stripped.dims.in = zeros(1,0);
    f_stripped.vars.in = cell(1,0);
    f_stripped.dom.in = zeros(0,2);
    f_stripped.depmat.in = false(1,0);
    f_stripped.chkval = false;

    f_init = initialize(f_stripped);
    if ~isequal(f_init.dims.in,d)
        error("Test_polyopvar_initialize failed: dims.in(i)>1 not correctly inferred (trial=%d,d=%d).",trial,d)
    end
    xval = rand_poly([d,1],xvars,2,4);
    Y_orig = apply_poly(f,{xval});
    Y_init = apply_poly(f_init,{xval});
    err = cleanpoly(Y_orig-Y_init,1e-10);
    if max(max(abs(err.C)))>tol
        error("Test_polyopvar_initialize failed: value changed after initialize for vector state (trial=%d).",trial)
    end

    % Case B (known limit, see NOTES): state only appears squared (no bare x),
    % so the scalar fallback in initialize.m is wrong for d>1; its final
    % consistency check (fdim_i~=Cdim_i) must catch this and error.
    if d>1
        g = x.*x;
        g_stripped = g;
        g_stripped.dims.in = zeros(1,0);
        g_stripped.vars.in = cell(1,0);
        g_stripped.dom.in = zeros(0,2);
        g_stripped.depmat.in = false(1,0);
        g_stripped.chkval = false;
        try
            initialize(g_stripped);
            error("Test_polyopvar_initialize failed: initialize() unexpectedly recovered dims.in(i)>1 for a state appearing only at degree>=2 -- the documented limitation (see header) may have been fixed; update this test and the report accordingly.")
        catch ME
            if contains(ME.message,'Test_polyopvar_initialize')
                rethrow(ME);
            end
        end
    end
end
disp('Test_polyopvar_initialize (vector-state) passed.')

% Validation checks: (1) non-tensopmat C, (2) dense Z, (3) wrong depmat.in,
% (4) duplicate vars.in.
x =polyopvar('x',{'s1'},[0,1]);

% (1) C not a 'tensopmat' object.
x_bad1 = x;             x_bad1.C = 5;           x_bad1.chkval = false;
try
    initialize(x_bad1);
    error("Test_polyopvar_initialize failed: non-tensopmat C was not rejected.")
catch ME
    if contains(ME.message,'Test_polyopvar_initialize'), rethrow(ME); end
    if ~contains(ME.message,'tensopmat')
        error("Test_polyopvar_initialize failed: non-tensopmat C was rejected, but not with an informative message (got: %s).",ME.message)
    end
end

% (2) Z coerced to sparse if supplied dense.
f2 = x.*x;
f2_bad = f2;            f2_bad.Z = full(f2.Z);  f2_bad.chkval = false;
f2_init = initialize(f2_bad);
if ~issparse(f2_init.Z)
    error("Test_polyopvar_initialize failed: dense Z was not coerced to sparse.")
end

% (3) Wrong depmat.in for a pure-power-only state (never appears bare).
f3 = x.*x;
f3_bad = f3;            f3_bad.depmat.in = ~f3.depmat.in;    f3_bad.chkval = false;
try
    initialize(f3_bad);
    error("Test_polyopvar_initialize failed: wrong depmat.in for a pure-power-only state was not rejected.")
catch ME
    if contains(ME.message,'Test_polyopvar_initialize'), rethrow(ME); end
end

% (4) Duplicate names in vars.in.
x4 = polyopvar('x',{'s1','s2'},[0,1;0,1]);
x4_bad = x4;
x4_bad.vars.in = {'s1','s1'};
x4_bad.dom.in = x4.dom.in([1,1],:);
x4_bad.chkval = false;
try
    initialize(x4_bad);
    error("Test_polyopvar_initialize failed: duplicate names in vars.in were not rejected.")
catch ME
    if contains(ME.message,'Test_polyopvar_initialize'), rethrow(ME); end
    if ~contains(ME.message,'unique')
        error("Test_polyopvar_initialize failed: duplicate vars.in names were rejected, but not with the specific, direct message (got: %s).",ME.message)
    end
end

disp('Test_polyopvar_initialize (validation-gap fixes) passed.')

function Test_polyopvar_initialize()
% Test_polyopvar_initialize.m -- one local function for each case:
%   case_basic              verify @polyopvar/initialize.m correctly infers
%                           missing fields (dims/vars/dom/depmat) from the
%                           coefficient operator, across varying numbers of
%                           spatial variables and
%   case_leftjustify        verify @polyopvar/initialize.m moves the
%                           nonzeros of each row of Z to the leading
%                           columns (order kept)
%   case_vars               verify @polyopvar/initialize.m when f.vars.in
%                           and f.C.vars.in differ
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_polyopvar_initialize_leftjustify,
%   Test_polyopvar_initialize_vars into this test, as cases;

cases = {@case_basic, @case_leftjustify, @case_vars};
nfail = 0;
for i=1:numel(cases)
    tstart = tic;
    try
        cases{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(cases{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(cases{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_polyopvar_initialize: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_polyopvar_initialize: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- verify @polyopvar/initialize.m
% correctly infers missing fields (dims/vars/dom/depmat) from the
% coefficient operator, across varying numbers of spatial variables and
% state dimensions, and correctly rejects an internally inconsistent
% object.
%
% Oracle: the reference object is built by the constructor (x.*x,
% x+x.*x) before its fields are stripped/corrupted; values are compared
% with apply_poly.
%

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
end

%% --------------------------------------------------------------------------
function case_leftjustify()
% case_leftjustify -- verify @polyopvar/initialize.m
% moves the nonzeros of each row of Z to the leading columns (order kept).
%
% Oracle: Z is scrambled by hand and checked structurally, and
% independently by evaluating the polynomial (apply_poly) before/after.
%

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
        error("case_leftjustify failed: Z not left-justified.")
    elseif ~issparse(g.Z)
        error("case_leftjustify failed: Z not sparse.")
    end
    Fg = apply_poly(g,{xv{1},yv{1}});
    err = cleanpoly(Fg-Fref,1e-10);
    if max(max(abs(err.C)))>tol*max(1,max(max(abs(Fref.C))))
        error("case_leftjustify failed: polynomial changed.")
    end

    % Already left-justified input is unchanged
    h = f;      h.chkval = false;       h = initialize(h);
    if ~isequal(h.Z,f.Z)
        error("case_leftjustify failed: justified Z altered.")
    end

    disp("case_leftjustify passed.")
end

%% --------------------------------------------------------------------------
function case_vars()
% case_vars -- verify @polyopvar/initialize.m when
% f.vars.in and f.C.vars.in differ.
%
%   Supported: variable of f that C lacks (no monomial may depend on it),
%   variable of C that f lacks (C may not act on it), and reordered
%   f.vars.in. Rejected: monomial depending on a variable C lacks; C acting
%   on a variable f lacks; wrong domains.
% Oracle: expected values are set by hand from the declarations; the
% polynomial is evaluated with apply_poly before/after.
%

    tol = 1e-8;
    pvar s1 s2 s3

    % Two states on different variables: x(s1), y(s1,s2), domains [0,1], [0,2]
    x = polyopvar('x',[s1,s2],[0,1;0,2]);
    y = polyopvar('y',[s1,s2],[0,1;0,2]);
    f = x + y + x.*y;
    if ~isequal(f.vars.in,{'s1','s2'}) || ~isequal(f.dom.in,[0,1;0,2]) || ...
            ~isequal(logical(f.depmat.in),true(2,2))
        error("case_vars failed: unexpected reference object.")
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
        error("case_vars failed: dom/depmat not in order of f.vars.in.")
    end

    % (2) f has an extra variable s3 that C lacks, no state depends on it
    g = f;      g.vars.in = {'s1','s2','s3'};
    g.dom.in = [0,1;0,2;0,3];       g.depmat.in = [f.depmat.in,false(2,1)];
    g.chkval = false;
    h = initialize(g);
    if ~h.chkval || ~isequal(h.vars.in,{'s1','s2','s3'}) || ~isequal(h.dom.in,[0,1;0,2;0,3])
        error("case_vars failed: extra variable of f not accepted.")
    end
    % ... and with depmat.in empty (inferred as no dependence on s3)
    g.depmat.in = false(0,3);       g.chkval = false;
    h = initialize(g);
    if ~isequal(logical(h.depmat.in),logical([1,1,0;1,1,0]))
        error("case_vars failed: depmat.in not inferred for extra variable.")
    end
    Fh = apply_poly(h,{xv,yv});
    err = cleanpoly(Fh-Fref,1e-10);
    if max(max(abs(err.C)))>tol*max(1,max(max(abs(Fref.C))))
        error("case_vars failed: polynomial changed.")
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
        error("case_vars failed: extra variable of C not accepted.")
    end
    g.vars.in = cell(1,0);  g.dom.in = zeros(0,2);  g.depmat.in = false(0,2);
    g.chkval = false;
    h = initialize(g);              % f.vars.in empty: copied from C
    if ~isequal(h.vars.in,{'s1','s2','s3'})
        error("case_vars failed: empty vars.in not copied from C.")
    end
    % 5b: C acts on s3
    C.depmat.in = [f.C.depmat.in,ones(size(f.C.depmat.in,1),1)];
    g = f;      g.C = C;            g.chkval = false;
    expect_error(@() initialize(g),'conflicts','C acting on variable missing in f');

    % (6) Wrong domain of a shared variable
    g = f;      g.dom.in = [0,1;0,5];       g.chkval = false;
    expect_error(@() initialize(g),'domains','wrong domain');

    disp("case_vars passed.")
end

%% ---- helpers --------------------------------------------------------------

function expect_error(fn,msgpart,label)
try
    fn();
catch ME
    if ~contains(ME.message,msgpart)
        error("case_vars failed (%s): rejected with unexpected message: %s",label,ME.message)
    end
    return
end
error("case_vars failed (%s): not rejected.",label)
end

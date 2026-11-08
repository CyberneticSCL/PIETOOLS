function results = test_convert_PIETOOLS_PDE_periodic()
%TEST_CONVERT_PIETOOLS_PDE_PERIODIC Exact, solver-free regressions of the
% conversion of 1D PDEs with periodic (and other) boundary conditions.
% Add PIETOOLS and the converters directory to the MATLAB path, then run
% results = test_convert_PIETOOLS_PDE_periodic. No SDP solver is required.
%
% The PIEs are 'copvar' operators on the components of the PIE state xf,
% see PIE.pde2pie. The tests evaluate the operators exactly on polynomials
% (pie_test_utils), and check the defining identities, instead of the number
% of auxiliary states, or the choice of their coordinates:
%   K*xf = 0,   T*xf = [auxiliary states; PDE states],
%   A*xf + B*inputs = [F*rhs; rhs],  with xf_b = C_b*D^{diff_b}x.
%
% DB, 11/07/2026: Rewritten for the 'copvar' PIEs of the converter. The tests
%   of the removed interface (PIE.F, PIE.misc.periodic, the argument for the
%   order of the components, the rejection of nearly singular BCs, which the
%   converter now converts) were replaced by the identities above.

tests = {@periodic_heat, @neumann_heat, @affine_kernel, @guaranteed_moment_fallback, ...
    @rank_zero_boundary_matrix, @mixed_components, @integral_boundaries, ...
    @shifted_domain_and_custom_variables, @repeated_boundary_terms, ...
    @inputs_and_outputs, @public_conversion_paths, ...
    @dirichlet_compatibility, @small_determinant_dirichlet, ...
    @nearly_singular_boundary, ...
    @first_order_transport, @coupled_ode_pde, ...
    @two_dimensional_dirichlet, @unsupported_boundary_forcing};
results = struct('name', cell(size(tests)), 'passed', false, 'message', '');
for k = 1:numel(tests)
    results(k).name = func2str(tests{k});
    try
        tests{k}();
        results(k).passed = true;
        fprintf('PASS %s\n', results(k).name);
    catch exception
        results(k).message = getReport(exception, 'extended', ...
            'hyperlinks', 'off');
        fprintf(2, 'FAIL %s\n%s\n', results(k).name, results(k).message);
    end
end
assert(all([results.passed]), 'PIETOOLS:PeriodicRegressionFailure', ...
    '%d of %d periodic conversion regressions failed.', ...
    nnz(~[results.passed]), numel(tests));
end

function periodic_heat()
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [1,0,-1,0; 0,1,0,-1]);
x = 3 + s^2*(1-s)^2;
[PIE, xf] = check_solution(PDE, x, 1);
check_constraint_detects_invalid_state(PIE, xf);
end

function neumann_heat()
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [0,1,0,0; 0,0,0,1]);
check_solution(PDE, 2 + 3*s^2*(1-s)^2, 1);
end

function affine_kernel()
% BCs x(0)=0 and x(1)-x_s(1)=0 leave the slope of an affine x free.
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [1,0,0,0; 0,0,1,-1]);
check_solution(PDE, 2*s + s^2*(1-s)^2, 1);
end

function guaranteed_moment_fallback()
% The heuristic moment 2*s is dependent on x(0)+(2/3)*x_s(0)=0.
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [1,2/3,0,0; 0,-1,0,1]);
check_solution(PDE, s-2/3+s^2*(1-s)^2, 1);
end

function rank_zero_boundary_matrix()
% Both BCs vanish on every affine function, so two moments are needed.
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [-1,-1,1,0; 0,-1,0,1]);
[PIE, xf] = check_solution(PDE, 2 + 3*s + s^2*(1-s)^2, 2);
check_constraint_detects_invalid_state(PIE, xf);
end

function mixed_components()
% One scalar state and one vector state, with affine nullities 1, 0, 2.
[s, sd] = variables();
n = 3;
B = zeros(2*n, 4*n);
B(1, [1, 2*n+1]) = [1, -1];
B(2, [n+1, 3*n+1]) = [1, -1];
B(3, 2) = 1;
B(4, 2*n+2) = 1;
B(5, [3, n+3, 2*n+3]) = [-1, -1, 1];
B(6, [n+3, 3*n+3]) = [-1, 1];
PDE = heat_pde([s, sd], [0, 1], [1, 2], B);
bubble = s^2*(1-s)^2;
x = [2 + bubble; s*(1-s); 4 + 2*s + 3*bubble];
check_solution(PDE, x, 3);
end

function integral_boundaries()
% Integral coefficients use the declared dummy variable. The first BC is
% integral_0^1 (theta-1/2)*x(theta) dtheta - x_s(0)/12 = 0.
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [0,0,0,0; 0,-1,0,1]);
PDE.BC{1}.term = {struct('x', 1, 'D', 0, 'loc', s, ...
    'I', {{[0,1]}}, 'C', sd-1/2), ...
    struct('x', 1, 'D', 1, 'loc', 0, 'C', -1/12)};
check_solution(PDE, 3 + 2*s + s^2*(1-s)^2, 2);
end

function shifted_domain_and_custom_variables()
s = pvar('periodic_space');
sd = pvar('periodic_dummy');
dom = [-2, 3];
r = s-dom(1);
L = diff(dom);
PDE = heat_pde([s, sd], dom, 1, [-1,-L,1,0; 0,-1,0,1]);
PIE = check_solution(PDE, 7 + 2*r + r^2*(L-r)^2, 2);
assert(isequal(double(PIE.space_pde.dom),dom), 'The domain must be preserved.');
end

function repeated_boundary_terms()
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [1,0,-1,0; 0,1,0,-1]);
% Splitting each coefficient across repeated terms must add contributions.
for k = 1:numel(PDE.BC)
    terms = PDE.BC{k}.term;
    for j = 1:numel(terms)
        terms{j}.C = terms{j}.C/3;
    end
    twice = terms;
    for j = 1:numel(twice)
        twice{j}.C = 2*twice{j}.C;
    end
    PDE.BC{k}.term = [terms, twice];
end
check_solution(PDE, 5 + 2*s^2*(1-s)^2, 1);
end

function inputs_and_outputs()
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [1,0,-1,0; 0,1,0,-1]);
PDE.x{1}.term{2} = struct('x', 1, 'D', 0, 'C', s+1);
PDE.x{1}.term{3} = struct('x', 1, 'D', 1, 'C', s^2);
PDE.u{1} = struct('vars', [], 'size', 1);
PDE.w{1} = struct('vars', [], 'size', 1);
PDE.x{1}.term{4} = struct('w', 1, 'C', 2-s);
PDE.x{1}.term{5} = struct('u', 1, 'C', 1+s);
PDE.z{1} = struct('vars', [], 'size', 1, 'term', ...
    {{struct('x', 1, 'loc', s, 'I', {{[0,1]}}, 'C', sd+1)}});
PDE.y{1} = struct('vars', [], 'size', 1, 'term', ...
    {{struct('x', 1, 'loc', 0, 'C', 1), ...
      struct('x', 1, 'loc', 1, 'D', 1, 'C', 2), ...
      struct('w', 1, 'C', 1), struct('u', 1, 'C', -1)}});
x = 3 + s^2*(1-s)^2;
w = 2;
u = -3;
rhs = diff(diff(x,s),s) + (s+1)*x + s^2*diff(x,s) + (2-s)*w + (1+s)*u;
PDE = initialize(PDE, true);
[PIE, xf] = check_solution(PDE, x, 1, rhs, {w}, {u});
[z, dzw, dzu] = deal(pie_test_utils('apply',PIE.Cz,xf), ...
    pie_test_utils('apply',PIE.Dzw,{w}), pie_test_utils('apply',PIE.Dzu,{u}));
assert_poly_equal(add_polys(z{1},dzw{1},dzu{1}), int((s+1)*x,s,0,1), 'integral output');
[y, dyw, dyu] = deal(pie_test_utils('apply',PIE.Cy,xf), ...
    pie_test_utils('apply',PIE.Dyw,{w}), pie_test_utils('apply',PIE.Dyu,{u}));
assert_poly_equal(add_polys(y{1},dyw{1},dyu{1}), subs(x,s,0)+2*subs(diff(x,s),s,1)+w-u, ...
    'boundary output with feedthrough');
end

function public_conversion_paths()
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [1,0,-1,0; 0,1,0,-1]);
assert(nargout('convert_PIETOOLS_PDE')==1, ...
    'The public PDE converter must return only the PIE structure.');
x = 2 + s^2*(1-s)^2;
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
PIE_public = convert(PDE);
PIE_initialized = initialize(PIE);
PIE_roundtrip = pie_struct(struct(PIE));
PIEs = {PIE, PIE_public, PIE_initialized, PIE_roundtrip};
PDE = initialize(PDE, true);
for k = 1:numel(PIEs)
    xf = pie_test_utils('xf',PIEs{k},PDE,x);
    Kx = pie_test_utils('apply',PIEs{k}.K,xf);
    assert(all(cellfun(@(c) poly_norm(c)<1e-8,Kx)), 'K*xf must vanish.');
    Tx = pie_test_utils('apply',PIEs{k}.T,xf);
    assert_poly_equal(Tx{end}, x, 'reconstruction of the PDE state');
end
PIE_top = convert_PIETOOLS_PDE(PDE, {'silent', 'Top'});
assert(~isempty(PIE_top.T), 'The flag Top must return T.');
display_text = evalc('display(PIE)');
assert(contains(display_text,'K'), 'The PIE display must show K.');
end

function dirichlet_compatibility()
[s, sd] = variables();
B = [1,0,0,0; 0,0,1,0];
PDE = heat_pde([s, sd], [0, 1], 1, B);
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
assert(isempty(PIE.K) || isempty(PIE.K.C), ...
    'Ordinary PDE conversion must leave K empty.');
x = s*(1-s)*(2+s);
PDE = initialize(PDE, true);
[aux, ~] = pie_test_utils('aux', PIE);
assert(isempty(aux), 'Dirichlet must not acquire auxiliary states.');
xf = pie_test_utils('xf',PIE,PDE,x);
Tx = pie_test_utils('apply',PIE.T,xf);
assert_poly_equal(Tx{1}, x, 'ordinary Dirichlet reconstruction');
% Row-scaled BCs give the same PIE
PDE_scaled = initialize(heat_pde([s, sd], [0, 1], 1, diag([1e-2, 1e2])*B), true);
PIE_scaled = convert_PIETOOLS_PDE(PDE_scaled, 'silent');
xf = pie_test_utils('xf',PIE_scaled,PDE_scaled,x);
Tx = pie_test_utils('apply',PIE_scaled.T,xf);
assert_poly_equal(Tx{1}, x, 'row-scaled Dirichlet reconstruction');
end

function small_determinant_dirichlet()
% Dirichlet BCs of a very small norm are not rejected: the PIE is the same.
[s, sd] = variables();
B = 1e-10*[1,0,0,0; 0,0,1,0];
PDE = initialize(heat_pde([s, sd], [0, 1], 1, B), true);
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
x = s*(1-s)*(2+s);
xf = pie_test_utils('xf',PIE,PDE,x);
Tx = pie_test_utils('apply',PIE.T,xf);
assert_poly_equal(Tx{end}, x, 'reconstruction for BCs of a small norm');
end

function nearly_singular_boundary()
% A nearly singular boundary matrix on a long interval is converted.
[s, sd] = variables();
B = [0,0,1e-10,0; 0,0,1,1e-7];
PDE = initialize(heat_pde([s, sd], [0, 1e6], 1, B), true);
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
assert(~isempty(PIE.T), 'The conversion must return a PIE.');
end

function first_order_transport()
[s, sd] = variables();
PDE = pde_struct();
PDE.vars = [s,sd];
PDE.dom = [0,1];
PDE.x{1} = struct('vars', [s,sd], 'dom', [0,1], 'size', 1, ...
    'diff', 1, 'term', {{struct('x',1,'D',1,'C',-1)}});
PDE.BC{1} = struct('term', {{struct('x',1,'loc',0)}});
PDE = initialize(PDE, true);
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
x = s+s^2;
xf = pie_test_utils('xf',PIE,PDE,x);
Tx = pie_test_utils('apply',PIE.T,xf);
assert_poly_equal(Tx{end}, x, 'transport reconstruction');
Ax = pie_test_utils('apply',PIE.A,xf);
assert_poly_equal(Ax{end}, -diff(x,s), 'transport dynamics');
end

function coupled_ode_pde()
[s, sd] = variables();
PDE = pde_struct();
PDE.vars = [s,sd];
PDE.dom = [0,1];
PDE.x{1} = struct('vars', [], 'size', 1, 'term', ...
    {{struct('x',1,'C',-1), struct('x',2,'D',1,'loc',0,'C',1)}});
PDE.x{2} = struct('vars', [s,sd], 'dom', [0,1], 'size', 1, ...
    'diff', 2, 'term', ...
    {{struct('x',2,'D',2,'C',1), struct('x',1,'C',s)}});
PDE.BC{1} = struct('term', {{struct('x',2,'loc',0)}});
PDE.BC{2} = struct('term', {{struct('x',2,'loc',1)}});
PDE = initialize(PDE, true);
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
q = 3;
xpde = s*(1-s)*(2+s);
x = [q; xpde];
rhs = [-q+subs(diff(xpde,s),s,0); diff(diff(xpde,s),s)+s*q];
check_identities(PIE, PDE, x, rhs, {}, {});
end

function two_dimensional_dirichlet()
s1 = pvar('s1'); s2 = pvar('s2');
t1 = pvar('s1_dum'); t2 = pvar('s2_dum');
PDE = pde_struct();
PDE.vars = [s1,t1;s2,t2];
PDE.dom = [0,1;0,1];
PDE.x{1} = struct('vars', PDE.vars, 'dom', PDE.dom, 'size', 1, ...
    'diff', [2,2], 'term', ...
    {{struct('x',1,'D',[2,0],'C',1), struct('x',1,'D',[0,2],'C',1)}});
locations = {[0,s2], [1,s2], [s1,0], [s1,1]};
for k = 1:4
    PDE.BC{k} = struct('term', {{struct('x',1,'loc',locations{k})}});
end
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
assert(isa(PIE.T,'copvar'), 'The 2D PIE must be a copvar.');
assert(isempty(PIE.K) || isempty(PIE.K.C), ...
    'Ordinary 2D PDE conversion must leave K empty.');
assert(isequal(size(PIE.T.C),[1,1]), 'Ordinary 2D Dirichlet has a single state.');
end

function unsupported_boundary_forcing()
[s, sd] = variables();
PDE = heat_pde([s,sd], [0,1], 1, [1,0,-1,0;0,1,0,-1]);
PDE.w{1} = struct('vars', [], 'size', 1);
PDE.BC{1}.term{end+1} = struct('w',1,'C',1);
try
    convert_PIETOOLS_PDE(PDE, 'silent');
catch exception
    assert(contains(exception.message,'inputs in the boundary conditions'), ...
        'Boundary forcing must fail with the documented unsupported-case error.');
    return
end
error('Boundary forcing with singular BCs must be rejected explicitly.');
end

function [PIE, xf] = check_solution(PDE, x, expected_aux, rhs, w, u)
% Convert, and check the identities of the PIE on the PDE state x
PDE = initialize(PDE, true);
assert_poly_equal(boundary_values(PDE,x), zeros(sum(PDE.BC_tab(:,2)),1), ...
    'test polynomial satisfies declared BCs');
if nargin<4 || isempty(rhs)
    rhs = x;
    for k = 1:numel(x)
        rhs(k) = diff(diff(x(k),PDE.vars(1,1)),PDE.vars(1,1));
    end
end
if nargin<5, w = {}; end
if nargin<6, u = {}; end
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
[aux, ~] = pie_test_utils('aux', PIE);
n_aux = sum(PIE.space.x.dims(aux));
assert(n_aux==expected_aux, ...
    'Unexpected number of auxiliary states: %d, expected %d.', n_aux, expected_aux);
xf = check_identities(PIE, PDE, x, rhs, w, u);
end

function xf = check_identities(PIE, PDE, x, rhs, w, u)
% K*xf=0, T*xf=[aux; PDE states], and A*xf+B*inputs=[F*rhs; rhs]
xf = pie_test_utils('xf', PIE, PDE, x);
p2p = PIE.pde2pie.x;
[aux, ~] = pie_test_utils('aux', PIE);
if ~isempty(PIE.K) && ~isempty(PIE.K.C)
    Kx = pie_test_utils('apply', PIE.K, xf);
    for r = 1:numel(Kx)
        assert(isempty(Kx{r}) || poly_norm(Kx{r})<1e-8, 'K*xf must vanish.');
    end
end
Tx = pie_test_utils('apply', PIE.T, xf);
Ax = pie_test_utils('apply', PIE.A, xf);
if ~isempty(w) && ~isempty(PIE.Bw), Bw = pie_test_utils('apply',PIE.Bw,w); else, Bw = []; end
if ~isempty(u) && ~isempty(PIE.Bu), Bu = pie_test_utils('apply',PIE.Bu,u); else, Bu = []; end
for r = 1:numel(xf)
    k = p2p.idx(r,1);
    rhs_k = pie_test_utils('pde_state', PDE, rhs, k);
    if ismember(r, aux)
        assert_poly_equal(Tx{r}, xf{r}, sprintf('auxiliary identity of block %d',r));
        expected = apply_F(PIE, r, rhs_k);
    else
        assert_poly_equal(Tx{r}, pie_test_utils('pde_state',PDE,x,k), ...
            sprintf('reconstruction of block %d',r));
        expected = rhs_k;
    end
    total = Ax{r};
    if ~isempty(Bw) && ~isempty(Bw{r}), total = total+Bw{r}; end
    if ~isempty(Bu) && ~isempty(Bu{r}), total = total+Bu{r}; end
    assert_poly_equal(total, expected, sprintf('dynamics of block %d',r));
end
end

function total = add_polys(varargin)
% Sum of polynomials, where [] is a zero operator output
total = [];
for k = 1:numel(varargin)
    if isempty(varargin{k}), continue, end
    if isempty(total), total = varargin{k}; else, total = total+varargin{k}; end
end
end

function q = apply_F(PIE, b, f)
% The map xf_b = C_b*(derivative), applied to f
Cb = PIE.pde2pie.x.C.C{b,b};
tmp = pie_test_utils('apply', copvar({Cb}), {f});
q = tmp{1};
end

function check_constraint_detects_invalid_state(PIE, xf)
% K must reject a distributed state that is not a derivative of a periodic one
[aux, main] = pie_test_utils('aux', PIE);
bad = xf;
for b = aux, bad{b} = polynomial(zeros(size(xf{b}))); end
for b = main, bad{b} = polynomial(ones(size(xf{b}))); end
Kx = pie_test_utils('apply', PIE.K, bad);
assert(any(cellfun(@(c) ~isempty(c) && poly_norm(c)>1e-7, Kx)), ...
    'K must reject an incompatible highest derivative.');
end

function PDE = heat_pde(vars, dom, sizes, B)
% B acts on [x(a); x_s(a); x(b); x_s(b)].
PDE = pde_struct();
PDE.vars = vars;
PDE.dom = dom;
n = sum(sizes);
offsets = [0,cumsum(sizes)];
for j = 1:numel(sizes)
    PDE.x{j} = struct('vars', vars, 'dom', dom, 'size', sizes(j), ...
        'diff', 2, 'term', {{struct('x', j, 'D', 2, 'C', eye(sizes(j)))}});
end
for i = 1:size(B,1)
    terms = {};
    for block = 1:4
        for j = 1:numel(sizes)
            C = B(i, (block-1)*n+(offsets(j)+1:offsets(j+1)));
            if any(C)
                terms{end+1} = struct('x', j, 'C', C, ... %#ok<AGROW>
                    'D', mod(block-1,2), 'loc', dom(1+(block>2)));
            end
        end
    end
    % A placeholder permits replacing a zero row with an integral BC.
    if isempty(terms)
        terms = {struct('x', 1, 'C', zeros(1,sizes(1)), 'loc', dom(1))};
    end
    PDE.BC{i} = struct('vars', [], 'size', 1, 'term', {terms});
end
end

function values = boundary_values(PDE, x)
s = PDE.vars(1,1);
sd = PDE.vars(1,2);
offsets = [0;cumsum(PDE.x_tab(:,2))];
values = polynomial(zeros(0,1));
for i = 1:numel(PDE.BC)
    value = polynomial(zeros(PDE.BC{i}.size,1));
    for j = 1:numel(PDE.BC{i}.term)
        term = PDE.BC{i}.term{j};
        state = x(offsets(term.x)+1:offsets(term.x+1));
        for order = 1:term.D
            state = diff(state,s);
        end
        if ~isempty(term.I{1})
            value = value + int(term.C*subs(state,s,sd),sd,term.I{1});
        else
            value = value + term.C*subs(state,s,term.loc);
        end
    end
    values = [values;value]; %#ok<AGROW>
end
end

function assert_poly_equal(actual, expected, label)
assert(isequal(size(actual),size(expected)), ...
    'Dimension mismatch in %s.', label);
error_size = poly_norm(polynomial(actual)-polynomial(expected));
assert(error_size <= 3e-8*(1+poly_norm(expected)), ...
    'Polynomial identity failed in %s (coefficient error %.3g).', label, error_size);
end

function value = poly_norm(p)
p = combine(polynomial(p));
coefficients = full(p.coefficient);
value = max([0;abs(coefficients(:))]);
end

function [s, sd] = variables()
s = pvar('s');
sd = pvar('theta');
end

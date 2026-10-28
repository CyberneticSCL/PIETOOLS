function results = test_convert_PIETOOLS_PDE_periodic()
% CODEX generated test
%TEST_CONVERT_PIETOOLS_PDE_PERIODIC Exact, solver-free PBC regressions.
% Add PIETOOLS and the updated converters directory to the MATLAB path, then
% run results = test_convert_PIETOOLS_PDE_periodic. No SDP solver or Symbolic
% Math Toolbox is required: the checks use PIETOOLS' polynomial class.
%
% The auxiliary coordinates need not use any particular nullspace basis.
% We test their defining identities using q = F*x and f = d_s^2*x instead.

% MMP, 10/28/2026: PIE.Fop, PIE.Kop and the struct fields 'Fop', 'Kop'
%                  read as F, K, the pie_struct property names since
%                  c578e16c; the old names raised "Unrecognized ... field
%                  'Fop'" in 13 of the 19 tests. The display check reads
%                  'F:', 'K:', as display.m prints them.

tests = {@periodic_heat, @neumann_heat, @affine_kernel, @guaranteed_moment_fallback, ...
    @rank_zero_boundary_matrix, @mixed_components, @integral_boundaries, ...
    @shifted_domain_and_custom_variables, @repeated_boundary_terms, ...
    @inputs_and_outputs, @public_conversion_paths, @reordered_source_metadata, ...
    @dirichlet_compatibility, @small_determinant_dirichlet_rejection, ...
    @nearly_singular_boundary_rejection, ...
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
[PIE, ~, K, q] = check_solution(PDE, x, 1);
check_constraint_detects_invalid_state(PIE, K, q, 1);
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
% This forces the guaranteed dual-moment construction from Lemma 5.
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [1,2/3,0,0; 0,-1,0,1]);
check_solution(PDE, s-2/3+s^2*(1-s)^2, 1);
end

function rank_zero_boundary_matrix()
% Both BCs vanish on every affine function, so rank(G)=0 and n_aux=2.
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [-1,-1,1,0; 0,-1,0,1]);
[PIE, ~, K, q] = check_solution(PDE, 2 + 3*s + s^2*(1-s)^2, 2);
check_constraint_detects_invalid_state(PIE, K, q, 1);
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
% Together with x_s(1)=x_s(0), this leaves both affine modes free.
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
assert(isequal(PIE.dom, dom));
assert(all(all(isequal(PIE.vars, [s, sd]))));
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
PDE = initialize(PDE, true);
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
% F = PIE.Fop;                                                              % MMP, 10/28/2026 (was)
F = PIE.F;                                                                  % MMP, 10/28/2026
% K = PIE.Kop;                                                              % MMP, 10/28/2026 (was)
K = PIE.K;                                                                  % MMP, 10/28/2026
[q, ~] = apply_pi(F, [], x);
f = diff(diff(x, s), s);
check_reconstruction(PIE, F, K, q, f, x);
w = 2;
u = -3;
rhs = f + (s+1)*x + s^2*diff(x,s) + (2-s)*w + (1+s)*u;
[a0, a1] = apply_pi(PIE.A, q, f);
[b10, b11] = apply_pi(PIE.B1, w, []);
[b20, b21] = apply_pi(PIE.B2, u, []);
[expected_qdot, ~] = apply_pi(F, [], rhs);
assert_poly_equal(a0+b10+b20, expected_qdot, 'input moment dynamics');
assert_poly_equal(a1+b11+b21, rhs, 'input PDE dynamics');
[z, ~] = apply_pi(PIE.C1, q, f);
[zw, ~] = apply_pi(PIE.D11, w, []);
[zu, ~] = apply_pi(PIE.D12, u, []);
assert_poly_equal(z+zw+zu, int((s+1)*x,s,0,1), 'integral output');
[y, ~] = apply_pi(PIE.C2, q, f);
[yw, ~] = apply_pi(PIE.D21, w, []);
[yu, ~] = apply_pi(PIE.D22, u, []);
assert_poly_equal(y+yw+yu, subs(x,s,0)+2*subs(diff(x,s),s,1)+w-u, ...
    'boundary output with feedthrough');
end

function public_conversion_paths()
[s, sd] = variables();
PDE = heat_pde([s, sd], [0, 1], 1, [1,0,-1,0; 0,1,0,-1]);
assert(nargout('convert_PIETOOLS_PDE')==1, ...
    'The public PDE converter must return only the PIE structure.');
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
% F = PIE.Fop;                                                              % MMP, 10/28/2026 (was)
F = PIE.F;                                                                  % MMP, 10/28/2026
% K = PIE.Kop;                                                              % MMP, 10/28/2026 (was)
K = PIE.K;                                                                  % MMP, 10/28/2026
PIE_initialized = initialize(PIE);
PIE_public = convert(PDE);
PIE_top = convert_PIETOOLS_PDE(PDE, {'silent', 'Top'});
% PIE_from_struct = pie_struct(struct('Fop', F, 'Kop', K));                 % MMP, 10/28/2026 (was)
PIE_from_struct = pie_struct(struct('F', F, 'K', K));                       % MMP, 10/28/2026
PIE_roundtrip = pie_struct(struct(PIE));
display_text = evalc('display(PIE)');
% assert(contains(display_text,'Fop:') && contains(display_text,'Kop:'), ...
%     'PIE display must show the auxiliary operators.');                    % MMP, 10/28/2026 (was)
assert(contains(display_text,'F:') && contains(display_text,'K:'), ...
    'PIE display must show the auxiliary operators.');                      % MMP, 10/28/2026
x = 2 + s^2*(1-s)^2;
[q, ~] = apply_pi(F, [], x);
f = diff(diff(x,s),s);
% check_reconstruction(PIE_initialized, PIE_initialized.Fop, ...
%     PIE_initialized.Kop, q, f, x);                                        % MMP, 10/28/2026 (was)
check_reconstruction(PIE_initialized, PIE_initialized.F, ...
    PIE_initialized.K, q, f, x);                                            % MMP, 10/28/2026
% check_reconstruction(PIE_public, PIE_public.Fop, PIE_public.Kop, q, f, x); % MMP, 10/28/2026 (was)
check_reconstruction(PIE_public, PIE_public.F, PIE_public.K, q, f, x);      % MMP, 10/28/2026
% check_reconstruction(PIE_top, PIE_top.Fop, PIE_top.Kop, q, f, x);         % MMP, 10/28/2026 (was)
check_reconstruction(PIE_top, PIE_top.F, PIE_top.K, q, f, x);               % MMP, 10/28/2026
% check_reconstruction(PIE_roundtrip, PIE_roundtrip.Fop, ...
%     PIE_roundtrip.Kop, q, f, x);                                          % MMP, 10/28/2026 (was)
check_reconstruction(PIE_roundtrip, PIE_roundtrip.F, ...
    PIE_roundtrip.K, q, f, x);                                              % MMP, 10/28/2026
% [copied_q, ~] = apply_pi(PIE_from_struct.Fop, [], x);                     % MMP, 10/28/2026 (was)
[copied_q, ~] = apply_pi(PIE_from_struct.F, [], x);                         % MMP, 10/28/2026
assert_poly_equal(copied_q, q, 'struct constructor preserves moment operator');
% [copied_k0, copied_k1] = apply_pi(PIE_from_struct.Kop, q, f);             % MMP, 10/28/2026 (was)
[copied_k0, copied_k1] = apply_pi(PIE_from_struct.K, q, f);                 % MMP, 10/28/2026
assert_poly_equal(copied_k0, zeros(K.dim(1,1),1), ...
    'struct constructor preserves finite compatibility constraint');
assert_poly_equal(copied_k1, zeros(K.dim(2,1),1), ...
    'struct constructor preserves distributed compatibility constraint');
% check_constraint_detects_invalid_state(PIE, PIE_from_struct.Kop, q, 1);   % MMP, 10/28/2026 (was)
check_constraint_detects_invalid_state(PIE, PIE_from_struct.K, q, 1);       % MMP, 10/28/2026
end

function dirichlet_compatibility()
[s, sd] = variables();
B = [1,0,0,0; 0,0,1,0];
PDE = heat_pde([s, sd], [0, 1], 1, B);
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
PIE_default = pie_struct();
% assert(isempty(PIE_default.Fop) && isempty(PIE_default.Kop), ...
%     'New PIE structures must default to empty auxiliary operators.');     % MMP, 10/28/2026 (was)
assert(isempty(PIE_default.F) && isempty(PIE_default.K), ...
    'New PIE structures must default to empty auxiliary operators.');       % MMP, 10/28/2026
% assert(isempty(PIE.Fop) && isempty(PIE.Kop), ...
%     'Ordinary PDE conversion must leave auxiliary operators empty.');     % MMP, 10/28/2026 (was)
assert(isempty(PIE.F) && isempty(PIE.K), ...
    'Ordinary PDE conversion must leave auxiliary operators empty.');       % MMP, 10/28/2026
x = s*(1-s)*(2+s);
f = diff(diff(x,s),s);
[t0, t1] = apply_pi(PIE.T, [], f);
assert(isempty(t0));
assert_poly_equal(t1, x, 'ordinary Dirichlet reconstruction');
assert(PIE.T.dim(1,1)==0, 'Dirichlet must not acquire auxiliary states.');
% Use scaling for which the original pinv retains both boundary directions.
PDE_scaled = heat_pde([s, sd], [0, 1], 1, diag([1e-2, 1e2])*B);
PIE_scaled = convert_PIETOOLS_PDE(PDE_scaled, 'silent');
[~, scaled_x] = apply_pi(PIE_scaled.T, [], f);
assert_poly_equal(scaled_x, x, 'row-scaled Dirichlet reconstruction');
assert(PIE_scaled.T.dim(1,1)==0);
end

function small_determinant_dirichlet_rejection()
% Preserve the original converter's rejection even though G is full rank.
% Its determinant is 1e-20 <= eps, and this must not select the PBC path.
[s, sd] = variables();
B = 1e-10*[1,0,0,0; 0,0,1,0];
PDE = heat_pde([s, sd], [0, 1], 1, B);
try
    convert_PIETOOLS_PDE(PDE, 'silent');
catch exception
    assert(contains(exception.message, ...
        'The PI operator defining the boundary conditions appears not to be invertible.'), ...
        'Small full-rank boundary matrices must keep the original invertibility error.');
    return
end
error('A full-rank boundary matrix with determinant <= eps must be rejected.');
end

function nearly_singular_boundary_rejection()
% Numerical rank alone would classify this matrix as singular, but its
% determinant is nonzero and safely below eps. A long interval separates
% the nonzero determinant from the numerical rank of the scaled matrix.
[s, sd] = variables();
B = [0,0,1e-10,0; 0,0,1,1e-7];
PDE = heat_pde([s, sd], [0, 1e6], 1, B);
try
    convert_PIETOOLS_PDE(PDE, 'silent');
catch exception
    assert(contains(exception.message, ...
        'The PI operator defining the boundary conditions appears not to be invertible.'), ...
        'Nearly singular full-rank boundary matrices must keep the original error.');
    return
end
error('A nonzero determinant <= eps must not select the periodic conversion.');
end

function reordered_source_metadata()
[s, sd] = variables();
B = zeros(4,8);
B(1,[1,5]) = [1,-1];
B(2,[3,7]) = [1,-1];
B(3,[2,4,6]) = [-1,-1,1];
B(4,[4,8]) = [-1,1];
PDE = heat_pde([s, sd], [0,1], [1,1], B);
order = struct('x', [2;1], 'y', [], 'z', [], 'u', [], 'w', []);
PIE = convert_PIETOOLS_PDE(PDE, order, {'silent'});
PIE_top = convert_PIETOOLS_PDE(PDE, order, {'silent','Top'});
assert(isequal(PIE.misc.periodic.source_components, ...
    PIE_top.misc.periodic.source_components), ...
    'Top and full conversion must preserve the same source component IDs.');
end

function first_order_transport()
[s, sd] = variables();
PDE = pde_struct();
PDE.vars = [s,sd];
PDE.dom = [0,1];
PDE.x{1} = struct('vars', [s,sd], 'dom', [0,1], 'size', 1, ...
    'diff', 1, 'term', {{struct('x',1,'D',1,'C',-1)}});
PDE.BC{1} = struct('term', {{struct('x',1,'loc',0)}});
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
x = s+s^2;
f = diff(x,s);
[t0,t1] = apply_pi(PIE.T, [], f);
assert(isempty(t0));
assert_poly_equal(t1, x, 'legacy transport reconstruction');
[a0,a1] = apply_pi(PIE.A, [], f);
assert(isempty(a0));
assert_poly_equal(a1, -f, 'legacy transport dynamics');
assert(isequal(PIE.T.dim,[0,0;1,1]));
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
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
x = s*(1-s)*(2+s);
q = 3;
f = diff(diff(x,s),s);
[t0,t1] = apply_pi(PIE.T, q, f);
assert_poly_equal(t0, q, 'legacy coupled ODE identity');
assert_poly_equal(t1, x, 'legacy coupled PDE reconstruction');
[a0,a1] = apply_pi(PIE.A, q, f);
assert_poly_equal(a0, -q+subs(diff(x,s),s,0), 'legacy coupled ODE dynamics');
assert_poly_equal(a1, f+s*q, 'legacy coupled PDE dynamics');
assert(isequal(PIE.T.dim,[1,1;1,1]));
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
assert(isa(PIE.T,'opvar2d'));
% assert(isempty(PIE.Fop) && isempty(PIE.Kop), ...
%     'Ordinary 2D PDE conversion must leave auxiliary operators empty.');  % MMP, 10/28/2026 (was)
assert(isempty(PIE.F) && isempty(PIE.K), ...
    'Ordinary 2D PDE conversion must leave auxiliary operators empty.');    % MMP, 10/28/2026
assert(isequal(PIE.T.dim,[0,0;0,0;0,0;1,1]));
x = s1*(1-s1)*s2*(1-s2);
f = diff(diff(diff(diff(x,s1),s1),s2),s2);
assert_poly_equal(apply_pure_2d(PIE.T,f), x, 'legacy 2D reconstruction');
rhs = diff(diff(x,s1),s1)+diff(diff(x,s2),s2);
assert_poly_equal(apply_pure_2d(PIE.A,f), rhs, 'legacy 2D heat dynamics');
end

function unsupported_boundary_forcing()
[s, sd] = variables();
PDE = heat_pde([s,sd], [0,1], 1, [1,0,-1,0;0,1,0,-1]);
PDE.w{1} = struct('vars', [], 'size', 1);
PDE.BC{1}.term{end+1} = struct('w',1,'C',1);
try
    convert_PIETOOLS_PDE(PDE, 'silent');
catch exception
    assert(strcmp(exception.identifier,'PIETOOLS:convert:SingularBCInputs'), ...
        'Boundary forcing must fail with the documented unsupported-case error.');
    return
end
error('Boundary forcing with singular BCs must be rejected explicitly.');
end

function [PIE, F, K, q] = check_solution(PDE, x, expected_aux)
PDE = initialize(PDE, true);
assert_poly_equal(boundary_values(PDE,x), zeros(sum(PDE.BC_tab(:,2)),1), ...
    'test polynomial satisfies declared BCs');
PIE = convert_PIETOOLS_PDE(PDE, 'silent');
% F = PIE.Fop;                                                              % MMP, 10/28/2026 (was)
F = PIE.F;                                                                  % MMP, 10/28/2026
% K = PIE.Kop;                                                              % MMP, 10/28/2026 (was)
K = PIE.K;                                                                  % MMP, 10/28/2026
assert(isfield(PIE.misc, 'periodic'), 'Missing periodic conversion metadata.');
assert(PIE.misc.periodic.n_aux==expected_aux, ...
    'Unexpected number of auxiliary moment states.');
assert(isequal(PIE.T.dim, [expected_aux,expected_aux;numel(x),numel(x)]), ...
    'Augmented PIE.T must be square.');
assert(sum(PIE.x_tab(PIE.x_tab(:,3)==0,2))==expected_aux, ...
    'The finite auxiliary states must appear in PIE.x_tab.');
assert(sum(PIE.x_tab(PIE.x_tab(:,3)==1,2))==numel(x));
[q, ~] = apply_pi(F, [], x);
assert(numel(q)==expected_aux);
s = PDE.vars(1,1);
f = diff(diff(x,s),s);
check_reconstruction(PIE, F, K, q, f, x);
[a0, a1] = apply_pi(PIE.A, q, f);
[moment_rhs, ~] = apply_pi(F, [], f);
assert_poly_equal(a0, moment_rhs, 'auxiliary heat dynamics');
assert_poly_equal(a1, f, 'distributed heat dynamics');
end

function check_reconstruction(PIE, F, K, q, f, x)
[r0, reconstructed] = apply_pi(PIE.misc.periodic.reconstruction, q, f);
assert(isempty(r0), 'PDE reconstruction must have no finite output.');
assert_poly_equal(reconstructed, x, 'PDE reconstruction');
s = PIE.vars(1,1);
assert_poly_equal(diff(diff(reconstructed,s),s), f, 'second derivative identity');
[recovered_q, ~] = apply_pi(F, [], reconstructed);
assert_poly_equal(recovered_q, q, 'moment coordinate identity');
[k0, k1] = apply_pi(K, q, f);
assert_poly_equal(k0, zeros(K.dim(1,1),1), 'finite compatibility constraint');
assert_poly_equal(k1, zeros(K.dim(2,1),1), 'distributed compatibility constraint');
[t0, t1] = apply_pi(PIE.T, q, f);
assert_poly_equal(t0, q, 'augmented finite identity');
assert_poly_equal(t1, x, 'augmented reconstruction');
end

function check_constraint_detects_invalid_state(PIE, K, q, f)
[k0, k1] = apply_pi(K, q, polynomial(f));
assert(poly_norm(k0)+poly_norm(k1)>1e-7, ...
    'K must reject an incompatible highest derivative.');
[~, invalid_x] = apply_pi(PIE.misc.periodic.reconstruction, q, polynomial(f));
assert(poly_norm(diff(diff(invalid_x,PIE.vars(1,1)),PIE.vars(1,1))-f)<1e-7);
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

function [finite, distributed] = apply_pi(P, q, f)
% Evaluate an opvar on an explicit polynomial, integrating its kernels.
assert(isa(P,'opvar'), 'Expected a one-dimensional opvar.');
q = polynomial(reshape(q,numel(q),1));
f = polynomial(reshape(f,numel(f),1));
assert(numel(q)==P.dim(1,2) && numel(f)==P.dim(2,2), ...
    'Polynomial test input dimensions do not match the operator.');
s = P.var1;
sd = P.var2;
finite = polynomial(zeros(P.dim(1,1),1));
distributed = polynomial(zeros(P.dim(2,1),1));
if P.dim(1,1)>0
    if P.dim(1,2)>0, finite = finite+P.P*q; end
    if P.dim(2,2)>0
        finite = finite+int(P.Q1*f,s,P.I(1),P.I(2));
    end
end
if P.dim(2,1)>0
    if P.dim(1,2)>0, distributed = distributed+P.Q2*q; end
    if P.dim(2,2)>0
        distributed = distributed+P.R.R0*f ...
            + int(P.R.R1*subs(f,s,sd),sd,P.I(1),s) ...
            + int(P.R.R2*subs(f,s,sd),sd,s,P.I(2));
    end
end
end

function result = apply_pure_2d(P, f)
% Evaluate the 2D-to-2D block when all other operator dimensions vanish.
assert(all(all(P.dim(1:3,:)==0)));
s = P.var1;
t = P.var2;
result = polynomial(zeros(P.dim(4,1),1));
for i = 1:3
    for j = 1:3
        term = f;
        if i>1, term = subs(term,s(1),t(1)); end
        if j>1, term = subs(term,s(2),t(2)); end
        term = P.R22{i,j}*term;
        if i==2, term = int(term,t(1),P.I(1,1),s(1)); end
        if i==3, term = int(term,t(1),s(1),P.I(1,2)); end
        if j==2, term = int(term,t(2),P.I(2,1),s(2)); end
        if j==3, term = int(term,t(2),s(2),P.I(2,2)); end
        result = result+term;
    end
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

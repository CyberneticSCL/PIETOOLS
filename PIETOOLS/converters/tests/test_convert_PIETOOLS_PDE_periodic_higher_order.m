function results = test_convert_PIETOOLS_PDE_periodic_higher_order()
%TEST_CONVERT_PIETOOLS_PDE_PERIODIC_HIGHER_ORDER Solver-free FTC regressions.
% Tests scalar, vector and mixed differentiability orders without assuming a
% particular basis for the auxiliary moments, by the identities of the PIE
% on the 'copvar' operators, see test_convert_PIETOOLS_PDE_periodic:
%   K*xf = 0,   T*xf = [auxiliary states; PDE states],
%   A*xf + B*inputs = [F*rhs; rhs].
% No SDP solver is required.
%
% DB, 11/07/2026: Rewritten for the 'copvar' PIEs of the converter, as for
%   test_convert_PIETOOLS_PDE_periodic (no PIE.F, PIE.misc.periodic, or
%   x_tab ordering; the identities above instead).

tests = {@periodic_first_order, @periodic_third_order, ...
    @periodic_fourth_order, @periodic_sixth_order, ...
    @guaranteed_moment_fourth_order, @guaranteed_moment_sixth_order, ...
    @rank_zero_fourth_order, @mixed_orders_and_vector_component, ...
    @mixed_zero_and_fourth_order, ...
    @shifted_domain_and_variables, @integral_third_order_boundaries, ...
    @higher_derivative_inputs_and_outputs, @public_conversion_paths, ...
    @nonsingular_fourth_order};
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
assert(all([results.passed]), 'PIETOOLS:HigherOrderPeriodicRegressionFailure', ...
    '%d of %d higher-order periodic regressions failed.', ...
    nnz(~[results.passed]), numel(tests));
end

function periodic_first_order()
check_scalar_periodic(1);
end

function periodic_third_order()
check_scalar_periodic(3);
end

function periodic_fourth_order()
check_scalar_periodic(4);
end

function periodic_sixth_order()
check_scalar_periodic(6);
end

function check_scalar_periodic(order)
[s, sd] = variables();
PDE = make_pde([s,sd], [0,1], order, 1, 'periodic');
PDE.x{1}.term{2} = struct('x',1,'D',0,'C',s+1);
if order>1
    PDE.x{1}.term{3} = struct('x',1,'D',order-1,'C',s);
end
x = 2+s^order*(1-s)^order;
[PIE,xf] = check_solution(PDE, x, 1);
check_rejects_constant_highest_derivative(PIE, xf);
end

function rank_zero_fourth_order()
% Taylor-matching endpoint BCs vanish on every polynomial of degree < 4.
[s, sd] = variables();
PDE = make_pde([s,sd], [0,1], 4, 1, 'taylor');
x = 2+3*s+4*s^2+5*s^3+s^4*(1-s)^4;
check_solution(PDE, x, 4);
end

function guaranteed_moment_fourth_order()
% A cubic moment cannot distinguish the free cubic kernel from the BCs.
check_guaranteed_moment(4,2/21);
end

function guaranteed_moment_sixth_order()
% The corresponding order-six fixture, with a fifth-degree kernel.
check_guaranteed_moment(6,1/220);
end

function check_guaranteed_moment(order,alpha)
[s, sd] = variables();
PDE = make_pde([s,sd],[0,1],order,1,'periodic');
PDE.BC{1}.term = {struct('x',1,'D',0,'loc',0,'C',1), ...
    struct('x',1,'D',order-1,'loc',0,'C',alpha)};
for k = 1:order-2
    PDE.BC{k+1}.term = {struct('x',1,'D',k,'loc',0,'C',1)};
end
x = s^(order-1)/factorial(order-1)-alpha+s^order*(1-s)^order;
check_solution(PDE,x,1);
end

function mixed_orders_and_vector_component()
% The first component is a two-vector of order 4, the second is scalar of
% order 2.
[s, sd] = variables();
PDE = make_pde([s,sd], [0,1], [4,2], [2,1], 'periodic');
PDE.x{1}.term{2} = struct('x',2,'D',1,'C',[1;2]);
PDE.x{2}.term{2} = struct('x',1,'D',3,'C',[2,-1]);
x = [2+s^4*(1-s)^4; 3+2*s^4*(1-s)^4; 4+3*s^2*(1-s)^2];
check_solution(PDE, x, 3);
end

function shifted_domain_and_variables()
s = pvar('higher_order_space');
sd = pvar('higher_order_dummy');
dom = [-1,1];
PDE = make_pde([s,sd], dom, 3, 1, 'periodic');
x = 3+(s-dom(1))^3*(dom(2)-s)^3;
PIE = check_solution(PDE, x, 1);
assert(isequal(double(PIE.space_pde.dom),dom), 'The domain must be preserved.');
end

function mixed_zero_and_fourth_order()
% A distributed state with no spatial differentiability is already a
% fundamental coordinate.
[s, sd] = variables();
PDE = make_pde([s,sd],[0,1],[4,0],[1,1],'periodic');
PDE.x{1}.term{2} = struct('x',2,'D',0,'C',s);
PDE.x{2}.term{2} = struct('x',1,'D',2,'C',1);
x = [2+s^4*(1-s)^4;3+s+s^2];
check_solution(PDE,x,1);
end

function integral_third_order_boundaries()
% All degree <= 2 kernel modes satisfy this integral BC and the remaining
% Taylor-matching rows. The integral is also zero for the symmetric bubble.
[s, sd] = variables();
PDE = make_pde([s,sd], [0,1], 3, 1, 'taylor');
PDE.BC{1}.term = {struct('x',1,'D',2,'loc',s, ...
    'I',{{[0,1]}},'C',sd-1/2)};
x = 2+3*s+4*s^2+s^3*(1-s)^3;
check_solution(PDE, x, 3);
end

function higher_derivative_inputs_and_outputs()
[s, sd] = variables();
PDE = make_pde([s,sd], [0,1], 4, 1, 'periodic');
PDE.x{1}.term{2} = struct('x',1,'D',3,'C',s^2);
PDE.x{1}.term{3} = struct('x',1,'D',1,'C',2-s);
PDE.x{1}.term{4} = struct('x',1,'D',0,'C',s+1);
PDE.u{1} = struct('vars',[],'size',1);
PDE.w{1} = struct('vars',[],'size',1);
PDE.x{1}.term{5} = struct('w',1,'C',2-s);
PDE.x{1}.term{6} = struct('u',1,'C',1+s);
PDE.z{1} = struct('vars',[],'size',1,'term', ...
    {{struct('x',1,'D',2,'loc',s,'I',{{[0,1]}},'C',sd+1)}});
PDE.y{1} = struct('vars',[],'size',1,'term', ...
    {{struct('x',1,'D',1,'loc',0,'C',1), ...
      struct('x',1,'D',3,'loc',1,'C',2), ...
      struct('w',1,'C',1),struct('u',1,'C',-1)}});
% This periodic polynomial has nonzero first and third endpoint derivatives.
x = 2+6*s^5-15*s^4+10*s^3-s+s^4*(1-s)^4;
w = 2; u = -3;
f = derivative(x,s,4);
rhs = f+s^2*derivative(x,s,3)+(2-s)*derivative(x,s,1) ...
    +(s+1)*x+(2-s)*w+(1+s)*u;
PDE = initialize(PDE,true);
[PIE,xf] = check_solution(PDE,x,1,rhs,{w},{u});
z = add_polys(first_of(pie_test_utils('apply',PIE.Cz,xf)), ...
    first_of(pie_test_utils('apply',PIE.Dzw,{w})), ...
    first_of(pie_test_utils('apply',PIE.Dzu,{u})));
assert_poly_equal(z,int((s+1)*derivative(x,s,2),s,0,1),'integral higher-derivative output');
y = add_polys(first_of(pie_test_utils('apply',PIE.Cy,xf)), ...
    first_of(pie_test_utils('apply',PIE.Dyw,{w})), ...
    first_of(pie_test_utils('apply',PIE.Dyu,{u})));
assert_poly_equal(y,subs(derivative(x,s,1),s,0) ...
    +2*subs(derivative(x,s,3),s,1)+w-u,'boundary higher-derivative output');
end

function public_conversion_paths()
[s, sd] = variables();
PDE = make_pde([s,sd],[0,1],4,1,'periodic');
x = 2+s^4*(1-s)^4;
PIE = convert_PIETOOLS_PDE(PDE,'silent');
PIE_initialized = initialize(PIE);
PDE = initialize(PDE,true);
rhs = state_rhs(PDE,x);
check_identities(PIE,PDE,x,rhs,{},{});
check_identities(PIE_initialized,PDE,x,rhs,{},{});
PIE_top = convert_PIETOOLS_PDE(PDE,{'silent','Top'});
assert(~isempty(PIE_top.T),'The flag Top must return T.');
end

function nonsingular_fourth_order()
% Clamped conditions retain the original nonsingular conversion path.
[s, sd] = variables();
PDE = initialize(make_pde([s,sd],[0,1],4,1,'clamped'),true);
PIE = convert_PIETOOLS_PDE(PDE,'silent');
assert(isempty(PIE.K) || isempty(PIE.K.C), ...
    'Nonsingular fourth-order conversion must not acquire auxiliary states.');
[aux,~] = pie_test_utils('aux',PIE);
assert(isempty(aux),'Nonsingular fourth-order conversion must not acquire auxiliary states.');
x = s^2*(1-s)^2;
check_identities(PIE,PDE,x,state_rhs(PDE,x),{},{});
end

function [PIE,xf] = check_solution(PDE,x,expected_aux,rhs,w,u)
PDE = initialize(PDE,true);
assert_poly_equal(boundary_values(PDE,x),zeros(sum(PDE.BC_tab(:,2)),1), ...
    'test polynomial satisfies all declared BCs');
if nargin<4 || isempty(rhs), rhs = state_rhs(PDE,x); end
if nargin<5, w = {}; end
if nargin<6, u = {}; end
PIE = convert_PIETOOLS_PDE(PDE,'silent');
[aux,~] = pie_test_utils('aux',PIE);
n_aux = sum(PIE.space.x.dims(aux));
assert(n_aux==expected_aux,'Unexpected auxiliary-state count: %d, expected %d.',n_aux,expected_aux);
xf = check_identities(PIE,PDE,x,rhs,w,u);
end

function check_rejects_constant_highest_derivative(PIE,xf)
[aux,main] = pie_test_utils('aux',PIE);
bad = xf;
for b = aux, bad{b} = polynomial(zeros(size(xf{b}))); end
for b = main, bad{b} = polynomial(ones(size(xf{b}))); end
Kx = pie_test_utils('apply',PIE.K,bad);
assert(any(cellfun(@(c) ~isempty(c) && poly_norm(c)>1e-7,Kx)), ...
    'Periodic BCs must reject a constant nonzero highest derivative.');
end

function v = first_of(c)
v = c{1};
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

function PDE = make_pde(vars,dom,orders,sizes,kind)
PDE = pde_struct(); PDE.vars = vars; PDE.dom = dom;
L = diff(dom);
bc = 0;
for j = 1:numel(orders)
    d = orders(j); n = sizes(j);
    PDE.x{j} = struct('vars',vars,'dom',dom,'size',n,'diff',d, ...
        'term',{{struct('x',j,'D',d,'C',(-1)^(j+1)*eye(n))}});
    for scalar = 1:n
        selector = zeros(1,n); selector(scalar) = 1;
        for k = 0:d-1
            bc = bc+1;
            switch kind
                case 'periodic'
                    terms = {struct('x',j,'D',k,'loc',dom(1),'C',selector), ...
                        struct('x',j,'D',k,'loc',dom(2),'C',-selector)};
                case 'taylor'
                    terms = {struct('x',j,'D',k,'loc',dom(2),'C',selector)};
                    for ell = k:d-1
                        terms{end+1} = struct('x',j,'D',ell,'loc',dom(1), ...
                            'C',-L^(ell-k)/factorial(ell-k)*selector); %#ok<AGROW>
                    end
                case 'clamped'
                    assert(d==4 && n==1);
                    terms = {struct('x',j,'D',mod(k,2), ...
                        'loc',dom(1+(k>=2)),'C',selector)};
                otherwise
                    error('Unknown fixture boundary type.');
            end
            PDE.BC{bc} = struct('vars',[],'size',1,'term',{terms});
        end
    end
end
end

function values = boundary_values(PDE,x)
s = PDE.vars(1,1); sd = PDE.vars(1,2);
offsets = [0;cumsum(PDE.x_tab(:,2))];
values = polynomial(zeros(0,1));
for i = 1:numel(PDE.BC)
    value = polynomial(zeros(PDE.BC{i}.size,1));
    for j = 1:numel(PDE.BC{i}.term)
        term = PDE.BC{i}.term{j};
        state = derivative(x(offsets(term.x)+1:offsets(term.x+1)),s,term.D);
        if ~isempty(term.I{1})
            value = value+int(term.C*subs(state,s,sd),sd,term.I{1});
        else
            value = value+term.C*subs(state,s,term.loc);
        end
    end
    values = [values;value]; %#ok<AGROW>
end
end

function rhs = state_rhs(PDE,x)
s = PDE.vars(1,1);
offsets = [0;cumsum(PDE.x_tab(:,2))];
rhs = polynomial(zeros(0,1));
for i = 1:numel(PDE.x)
    value = polynomial(zeros(PDE.x{i}.size,1));
    for j = 1:numel(PDE.x{i}.term)
        term = PDE.x{i}.term{j};
        assert(isfield(term,'x'),'Use explicit dynamics for input fixtures.');
        state = x(offsets(term.x)+1:offsets(term.x+1));
        value = value+term.C*derivative(state,s,term.D);
    end
    rhs = [rhs;value]; %#ok<AGROW>
end
end

function p = derivative(p,s,order)
for k = 1:order
    p = diff(p,s);
end
end

function assert_poly_equal(actual,expected,label)
assert(isequal(size(actual),size(expected)),'Dimension mismatch in %s.',label);
error_size = poly_norm(polynomial(actual)-polynomial(expected));
assert(error_size<=3e-8*(1+poly_norm(expected)), ...
    'Polynomial identity failed in %s (coefficient error %.3g).',label,error_size);
end

function value = poly_norm(p)
p = combine(polynomial(p));
coefficients = full(p.coefficient);
value = max([0;abs(coefficients(:))]);
end

function [s,sd] = variables()
s = pvar('s'); sd = pvar('theta');
end

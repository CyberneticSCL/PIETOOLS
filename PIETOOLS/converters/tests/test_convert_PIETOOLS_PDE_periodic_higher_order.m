function results = test_convert_PIETOOLS_PDE_periodic_higher_order()
%TEST_CONVERT_PIETOOLS_PDE_PERIODIC_HIGHER_ORDER Solver-free FTC regressions.
% Tests scalar, vector and mixed differentiability orders without assuming a
% particular basis for the auxiliary moments. No SDP solver is required.

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
PIE = check_solution(PDE, x, 1);
check_inverse_identities(PIE, order);
[k0,k1] = apply_pi(PIE.Kop, 0, 1);
assert(poly_norm(k0)+poly_norm(k1)>1e-7, ...
    'Periodic BCs must reject a constant nonzero highest derivative.');
end

function rank_zero_fourth_order()
% Taylor-matching endpoint BCs vanish on every polynomial of degree < 4.
[s, sd] = variables();
PDE = make_pde([s,sd], [0,1], 4, 1, 'taylor');
x = 2+3*s+4*s^2+5*s^3+s^4*(1-s)^4;
PIE = check_solution(PDE, x, 4);
assert(PIE.misc.periodic.boundary_rank==0);
check_inverse_identities(PIE, 4);
end

function guaranteed_moment_fourth_order()
% A cubic moment cannot distinguish the free cubic kernel from the BCs.
% Conversion must use independent dual moments instead.
check_guaranteed_moment(4,2/21);
end

function guaranteed_moment_sixth_order()
% The corresponding order-six fixture exercises the same fallback with
% a fifth-degree kernel and a higher-degree FTC integral.
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
PIE = check_solution(PDE,x,1);
check_inverse_identities(PIE,order);
end

function mixed_orders_and_vector_component()
% The original first component is a two-vector of order 4, while the
% second is scalar of order 2. Conversion reorders the latter first.
[s, sd] = variables();
PDE = make_pde([s,sd], [0,1], [4,2], [2,1], 'periodic');
PDE.x{1}.term{2} = struct('x',2,'D',1,'C',[1;2]);
PDE.x{2}.term{2} = struct('x',1,'D',3,'C',[2,-1]);
x = [2+s^4*(1-s)^4; 3+2*s^4*(1-s)^4; 4+3*s^2*(1-s)^2];
PIE = check_solution(PDE, x, 3);
physical = PIE.x_tab(PIE.x_tab(:,3)==1,:);
assert(isequal(physical(:,1),[2;1]), ...
    'Mixed differentiability orders must retain the original component IDs.');
assert(isequal(physical(:,2:4),[1,1,2;2,1,4]), ...
    'Physical metadata must preserve vector sizes and spatial orders.');
assert(isequal(sort(PIE.misc.periodic.source_components),[1;1;2]), ...
    'Every periodic scalar coordinate needs a correctly attributed moment.');
check_inverse_identities(PIE, [2;4;4]);
end

function shifted_domain_and_variables()
s = pvar('higher_order_space');
sd = pvar('higher_order_dummy');
dom = [-1,1];
PDE = make_pde([s,sd], dom, 3, 1, 'periodic');
x = 3+(s-dom(1))^3*(dom(2)-s)^3;
PIE = check_solution(PDE, x, 1);
assert(isequal(PIE.dom,dom));
assert(all(all(isequal(PIE.vars,[s,sd]))));
check_inverse_identities(PIE, 3);
end

function mixed_zero_and_fourth_order()
% A distributed state with no spatial differentiability is already a
% fundamental coordinate. Its reconstruction needs the R0 identity block.
[s, sd] = variables();
PDE = make_pde([s,sd],[0,1],[4,0],[1,1],'periodic');
PDE.x{1}.term{2} = struct('x',2,'D',0,'C',s);
PDE.x{2}.term{2} = struct('x',1,'D',2,'C',1);
x = [2+s^4*(1-s)^4;3+s+s^2];
PIE = check_solution(PDE,x,1);
physical = PIE.x_tab(PIE.x_tab(:,3)==1,:);
assert(isequal(physical(:,1),[2;1]));
assert(isequal(physical(:,4),[0;4]));
assert(isequal(PIE.misc.periodic.source_components,1));
check_inverse_identities(PIE,[0;4]);
end

function integral_third_order_boundaries()
% All degree <= 2 kernel modes satisfy this integral BC and the remaining
% Taylor-matching rows. The integral is also zero for the symmetric bubble.
[s, sd] = variables();
PDE = make_pde([s,sd], [0,1], 3, 1, 'taylor');
PDE.BC{1}.term = {struct('x',1,'D',2,'loc',s, ...
    'I',{{[0,1]}},'C',sd-1/2)};
x = 2+3*s+4*s^2+s^3*(1-s)^3;
PIE = check_solution(PDE, x, 3);
check_inverse_identities(PIE, 3);
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
PDE = initialize(PDE,true);
assert_poly_equal(boundary_values(PDE,x),zeros(4,1),'input/output fixture BCs');
PIE = convert_PIETOOLS_PDE(PDE,'silent');
[q,~] = apply_pi(PIE.Fop,[],x);
f = derivative(x,s,4);
check_reconstruction(PIE,q,f,x,4);
w = 2; u = -3;
rhs = f+s^2*derivative(x,s,3)+(2-s)*derivative(x,s,1) ...
    +(s+1)*x+(2-s)*w+(1+s)*u;
[a0,a1] = apply_pi(PIE.A,q,f);
[b10,b11] = apply_pi(PIE.B1,w,[]);
[b20,b21] = apply_pi(PIE.B2,u,[]);
[expected_qdot,~] = apply_pi(PIE.Fop,[],rhs);
assert_poly_equal(a0+b10+b20,expected_qdot,'higher-order input moment dynamics');
assert_poly_equal(a1+b11+b21,rhs,'higher-order input PDE dynamics');
[z,~] = apply_pi(PIE.C1,q,f);
[zw,~] = apply_pi(PIE.D11,w,[]);
[zu,~] = apply_pi(PIE.D12,u,[]);
assert_poly_equal(z+zw+zu,int((s+1)*derivative(x,s,2),s,0,1), ...
    'integral higher-derivative output');
[y,~] = apply_pi(PIE.C2,q,f);
[yw,~] = apply_pi(PIE.D21,w,[]);
[yu,~] = apply_pi(PIE.D22,u,[]);
assert_poly_equal(y+yw+yu,subs(derivative(x,s,1),s,0) ...
    +2*subs(derivative(x,s,3),s,1)+w-u,'boundary higher-derivative output');
end

function public_conversion_paths()
[s, sd] = variables();
PDE = make_pde([s,sd],[0,1],4,1,'periodic');
x = 2+s^4*(1-s)^4;
PIE = convert_PIETOOLS_PDE(PDE,'silent');
PIE_top = convert_PIETOOLS_PDE(PDE,{'silent','Top'});
PIE_initialized = initialize(PIE);
[q,~] = apply_pi(PIE.Fop,[],x);
f = derivative(x,s,4);
check_reconstruction(PIE_top,q,f,x,4);
check_reconstruction(PIE_initialized,q,f,x,4);
assert(isequal(PIE.misc.periodic.source_components, ...
    PIE_top.misc.periodic.source_components));
end

function nonsingular_fourth_order()
% Clamped conditions retain the original nonsingular conversion path.
[s, sd] = variables();
PDE = make_pde([s,sd],[0,1],4,1,'clamped');
PIE = convert_PIETOOLS_PDE(PDE,'silent');
assert(isempty(PIE.Fop) && isempty(PIE.Kop), ...
    'Nonsingular fourth-order conversion must not acquire auxiliary states.');
x = s^2*(1-s)^2;
f = derivative(x,s,4);
[t0,t1] = apply_pi(PIE.T,[],f);
assert(isempty(t0));
assert_poly_equal(t1,x,'nonsingular fourth-order reconstruction');
[a0,a1] = apply_pi(PIE.A,[],f);
assert(isempty(a0));
assert_poly_equal(a1,f,'nonsingular fourth-order dynamics');
end

function PIE = check_solution(PDE,x,expected_aux)
PDE = initialize(PDE,true);
assert_poly_equal(boundary_values(PDE,x),zeros(sum(PDE.BC_tab(:,2)),1), ...
    'test polynomial satisfies all declared BCs');
rhs = state_rhs(PDE,x);
PIE = convert_PIETOOLS_PDE(PDE,'silent');
assert(isfield(PIE.misc,'periodic'));
assert(PIE.misc.periodic.n_aux==expected_aux,'Unexpected auxiliary-state count.');
assert(isequal(PIE.T.dim,[expected_aux,expected_aux;numel(x),numel(x)]), ...
    'Augmented reconstruction must be square.');
finite_rows = PIE.x_tab(:,3)==0;
assert(sum(PIE.x_tab(finite_rows,2))==expected_aux);
assert(isequal(PIE.x_tab(finite_rows,1),PIE.misc.periodic.source_components));
[x,orders] = physical_order(PDE,PIE,x);
[rhs,~] = physical_order(PDE,PIE,rhs);
[q,~] = apply_pi(PIE.Fop,[],x);
f = highest_derivatives(x,PIE.vars(1,1),orders);
check_reconstruction(PIE,q,f,x,orders);
[a0,a1] = apply_pi(PIE.A,q,f);
[expected_qdot,~] = apply_pi(PIE.Fop,[],rhs);
assert_poly_equal(a0,expected_qdot,'auxiliary higher-order dynamics');
assert_poly_equal(a1,rhs,'distributed higher-order dynamics');
end

function check_reconstruction(PIE,q,f,x,orders)
[r0,reconstructed] = apply_pi(PIE.misc.periodic.reconstruction,q,f);
assert(isempty(r0));
assert_poly_equal(reconstructed,x,'physical reconstruction');
assert_poly_equal(highest_derivatives(reconstructed,PIE.vars(1,1),orders), ...
    f,'componentwise highest-derivative identity');
[recovered_q,~] = apply_pi(PIE.Fop,[],reconstructed);
assert_poly_equal(recovered_q,q,'moment identity');
[k0,k1] = apply_pi(PIE.Kop,q,f);
assert_poly_equal(k0,zeros(PIE.Kop.dim(1,1),1),'finite BC compatibility');
assert_poly_equal(k1,zeros(PIE.Kop.dim(2,1),1),'distributed BC compatibility');
[t0,t1] = apply_pi(PIE.T,q,f);
assert_poly_equal(t0,q,'augmented finite identity');
assert_poly_equal(t1,x,'augmented distributed reconstruction');
end

function check_inverse_identities(PIE,orders)
% These two inverse identities must hold even when f violates the BCs.
s = PIE.vars(1,1);
q = (1:PIE.misc.periodic.n_aux)';
n = PIE.T.dim(2,2);
f = polynomial(ones(n,1))+s*(1:n)';
[~,x] = apply_pi(PIE.misc.periodic.reconstruction,q,f);
assert_poly_equal(highest_derivatives(x,s,orders),f, ...
    'highest derivative of arbitrary reconstruction');
[recovered_q,~] = apply_pi(PIE.Fop,[],x);
assert_poly_equal(recovered_q,q,'moments of arbitrary reconstruction');
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

function [ordered,orders] = physical_order(PDE,PIE,value)
offsets = [0;cumsum(PDE.x_tab(:,2))];
rows = PIE.x_tab(PIE.x_tab(:,3)==1,:);
ordered = polynomial(zeros(0,1)); orders = zeros(0,1);
for k = 1:size(rows,1)
    j = rows(k,1);
    ordered = [ordered;value(offsets(j)+1:offsets(j+1))]; %#ok<AGROW>
    orders = [orders;repmat(rows(k,4),rows(k,2),1)]; %#ok<AGROW>
end
end

function f = highest_derivatives(x,s,orders)
orders = orders(:);
assert(numel(x)==numel(orders));
f = polynomial(zeros(numel(x),1));
for k = 1:numel(x)
    f(k) = derivative(x(k),s,orders(k));
end
end

function p = derivative(p,s,order)
for k = 1:order
    p = diff(p,s);
end
end

function [finite,distributed] = apply_pi(P,q,f)
% Evaluate opvar kernels by polynomial integration.
assert(isa(P,'opvar'));
q = polynomial(reshape(q,numel(q),1));
f = polynomial(reshape(f,numel(f),1));
assert(numel(q)==P.dim(1,2) && numel(f)==P.dim(2,2));
s = P.var1; sd = P.var2;
finite = polynomial(zeros(P.dim(1,1),1));
distributed = polynomial(zeros(P.dim(2,1),1));
if P.dim(1,1)>0
    if P.dim(1,2)>0, finite = finite+P.P*q; end
    if P.dim(2,2)>0, finite = finite+int(P.Q1*f,s,P.I(1),P.I(2)); end
end
if P.dim(2,1)>0
    if P.dim(1,2)>0, distributed = distributed+P.Q2*q; end
    if P.dim(2,2)>0
        distributed = distributed+P.R.R0*f ...
            +int(P.R.R1*subs(f,s,sd),sd,P.I(1),s) ...
            +int(P.R.R2*subs(f,s,sd),sd,s,P.I(2));
    end
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

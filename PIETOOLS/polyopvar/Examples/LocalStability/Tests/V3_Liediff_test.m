%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - V3_Liediff_test.m
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% CR, 10/06/2026: Initial coding.
% CR, 10/06/2026: Added d=1 and d=2 analytical Lie-derivative
%                  regressions based on Lem. 13.
% CR, 10/06/2026: Added a numerical central finite-difference regression.
% CR, 10/06/2026: Added block-wise adjoint-symmetry constraint tests.

% Test the basic construction, distributed polynomial degrees, and the
% analytical d=1 and d=2 Lie derivatives returned by V3_Liediff.m.

clearvars; clear stateNameGenerator; close all; clc;

% Add PIETOOLS to the path when this test is run directly from its Tests
% directory or by using MATLAB's run command.
test_dir      = fileparts(mfilename('fullpath'));
local_dir     = fileparts(test_dir);
pietools_root = fileparts(fileparts(fileparts(local_dir)));
addpath(genpath(pietools_root));


%% 1. Construct a scalar polynomial PIE for the tests.

% Fisher equation with linear and quadratic terms. Consequently, the PIE
% right-hand side has distributed monomial degrees l in {1,2}.
pvar s t
dom = [0,1];
u   = pde_var(s,dom);
alp =  5;
bet = -1;
PDE = [diff(u,t)==diff(u,s,2) + alp*u - bet*u^2;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];

PIE = convert(PDE);
f   = PIE.f;
x   = f.vartab;

% Test low-order cases that exercise both the d=1 multiplier basis and the
% d>1 tensor-product construction. Each row is [d,opdeg]. Positive opdeg
% values are used because composition of the current degree-zero d=1 basis
% with Top fails in the underlying nopvar/mtimes implementation.
test_cases = [1,1;
              1,2;
              2,1];

f_degs = unique(sum(f.degmat,2));


%% 2. Run the construction and degree checks.

fprintf('\n --- Running V3_Liediff construction and degree tests ---\n');

for test_idx = 1:size(test_cases,1)
    d     = test_cases(test_idx,1);
    opdeg = test_cases(test_idx,2);

    % Use a fresh program in every case so the Gram variables and symmetry
    % constraints from one test do not affect any subsequent test.
    clear stateNameGenerator
    prog = piesos_program(x);
    [prog,V3,dV3,Pcell] = V3_Liediff(prog,PIE,d,opdeg);

    
    %% Basic construction checks.
    assert(isa(V3,'polyopvar'), ...
        'V3_Liediff_test:InvalidV3Class', ...
        'V3 must be returned as a polyopvar object.');
    assert(isa(dV3,'polyopvar'), ...
        'V3_Liediff_test:InvaliddV3Class', ...
        'dV3 must be returned as a polyopvar object.');
    assert(isequal(size(V3),[1,1]), ...
        'V3_Liediff_test:InvalidV3Size', ...
        'V3 must be a scalar FDP.');
    assert(isequal(size(dV3),[1,1]), ...
        'V3_Liediff_test:InvaliddV3Size', ...
        'dV3 must be a scalar FDP.');
    assert(isequal(V3.varname,f.varname), ...
        'V3_Liediff_test:InvalidV3Variables', ...
        'V3 must be expressed in the fundamental PIE state.');
    assert(isequal(dV3.varname,f.varname), ...
        'V3_Liediff_test:InvaliddV3Variables', ...
        'dV3 must be expressed in the fundamental PIE state.');

    % The independent dummy states used to impose adjoint symmetry should
    % occur only in the program constraints, not in either returned FDP.
    assert(~any(contains(V3.varname,{'_left','_right'})), ...
        'V3_Liediff_test:DummyVariableInV3', ...
        'A symmetry-constraint dummy state appears in V3.');
    assert(~any(contains(dV3.varname,{'_left','_right'})), ...
        'V3_Liediff_test:DummyVariableIndV3', ...
        'A symmetry-constraint dummy state appears in dV3.');

    %% Distributed-degree checks.

    % V3 contains block (i,j) terms of distributed degree i+j. A term of
    % PIE degree l substituted into one differentiated factor produces a
    % dV3 term of degree i+j-1+l.
    V3_degs_expected  = [];
    dV3_degs_expected = [];
    for i = 1:d
        for j = i:d
            V3_degs_expected(end+1,1) = i+j; %#ok<SAGROW>
            for l_idx = 1:numel(f_degs)
                dV3_degs_expected(end+1,1) = ...
                    i+j-1+f_degs(l_idx); %#ok<SAGROW>
            end
        end
    end

    V3_degs_expected  = unique(V3_degs_expected);
    dV3_degs_expected = unique(dV3_degs_expected);
    V3_degs_actual    = unique(sum(V3.degmat,2));
    dV3_degs_actual   = unique(sum(dV3.degmat,2));

    assert(isequal(V3_degs_actual(:),V3_degs_expected(:)), ...
        'V3_Liediff_test:IncorrectV3Degrees', ...
        'The distributed degrees of V3 do not match i+j.');
    assert(isequal(dV3_degs_actual(:),dV3_degs_expected(:)), ...
        'V3_Liediff_test:IncorrectdV3Degrees', ...
        'The distributed degrees of dV3 do not match i+j-1+l.');

    %% Analytical Lie-derivative regressions for d=1 and d=2.

    % Reset the generated dummy-variable names, then independently rebuild
    % the U*x and U*T*x factors used in Lem. 13. This lets the analytical
    % expressions below use the returned Gram blocks Pcell without relying
    % on the dV3 assembled internally by V3_Liediff.
    clear stateNameGenerator
    [Zx,ZTx] = build_regression_factors(PIE,d,opdeg);

    % Construct the tensor powers needed by both the analytical derivative
    % and the independent symmetry-constraint regression.
    Zs1 = cell(d,1);
    Zs2 = cell(d,1);
    Zs1{1} = Zx;
    Zs2{1} = ZTx;
    for degree_idx = 2:d
        Zs1{degree_idx} = DMB(Zs1{degree_idx-1},Zx);
        Zs2{degree_idx} = DMB(Zs2{degree_idx-1},ZTx);
    end

    if d == 1
        % Lem. 13 with q=1:
        %
        % dV3 = 2*sum_l <U*x,Q_11*U*C_l*x^(otimes l)>.
        R11 = innerprod_v2(Zx,Zx,Pcell{1,1},'skip_combine');
        dV3_analytical = 2*substitute_pie_rhs(R11,2,f);
    elseif d == 2
        % Write the q=2 sum in Lem. 13 explicitly. The six expressions
        % below correspond respectively to (i,j,k) equal to
        % (1,1,1), (1,2,1), (1,2,2), (2,1,1), (2,2,1), and (2,2,2),
        % with the two (2,2,k) contributions appearing as the last two
        % summands.
        Zx2  = Zs1{2};
        ZxZTx = DMB(Zx,ZTx);
        ZTxZx = DMB(ZTx,Zx);

        R11  = innerprod_v2(Zx,Zx,Pcell{1,1},'skip_combine');
        R121 = innerprod_v2(Zx,ZxZTx,Pcell{1,2},'skip_combine');
        R122 = innerprod_v2(Zx,ZTxZx,Pcell{1,2},'skip_combine');
        R211 = innerprod_v2(Zx2,Zx,Pcell{2,1},'skip_combine');
        R221 = innerprod_v2(Zx2,ZxZTx,Pcell{2,2},'skip_combine');
        R222 = innerprod_v2(Zx2,ZTxZx,Pcell{2,2},'skip_combine');

        dV3_analytical = 2*( ...
            substitute_pie_rhs(R11, 2,f) + ...
            substitute_pie_rhs(R121,2,f) + ...
            substitute_pie_rhs(R122,3,f) + ...
            substitute_pie_rhs(R211,3,f) + ...
            substitute_pie_rhs(R221,3,f) + ...
            substitute_pie_rhs(R222,4,f));
    else
        error('V3_Liediff_test:MissingAnalyticalRegression', ...
            'No analytical Lie-derivative regression is defined for d=%d.',d);
    end

    verify_fdp_equality(dV3,dV3_analytical,d);

    %% Adjoint-symmetry constraint regression.

    % Reconstruct Qhat_ij^*-Qhat_ji using independent left and right
    % distributed states. For each i<=j, verify that every coefficient
    % equation in this operator identity lies in the row space of the
    % equality constraints returned in prog. The diagonal cases test
    % self-adjointness; the strict upper triangle tests the paired blocks.
    state_name = f.varname{1};
    left_name  = [state_name,'_left'];
    right_name = [state_name,'_right'];
    Zs1_left   = Zs1;
    Zs1_right  = Zs1;
    Zs2_left   = Zs2;
    Zs2_right  = Zs2;
    symmetry_error = 0;

    for i = 1:d
        Zs1_left{i}.varname  = {left_name};
        Zs1_right{i}.varname = {right_name};
        Zs2_left{i}.varname  = {left_name};
        Zs2_right{i}.varname = {right_name};

        for j = i:d
            Qhat_ij_adjoint = innerprod_v2( ...
                Zs2_left{j},Zs1_right{i},Pcell{i,j}');
            Qhat_ji = innerprod_v2( ...
                Zs1_left{j},Zs2_right{i},Pcell{j,i});
            symmetry_residual = Qhat_ij_adjoint-Qhat_ji;
            block_error = verify_symmetry_is_enforced( ...
                prog,x,symmetry_residual,i,j);
            symmetry_error = max(symmetry_error,block_error);
        end
    end

    %% Numerical finite-difference regression.

    % Select a deterministic, nonzero point in the nullspace of the
    % adjoint-symmetry equality constraints added by V3_Liediff. The
    % identity in Lem. 13 is only valid subject to these constraints.
    [dvar_names,dvar_values] = symmetry_feasible_point(prog);

    % For the Fisher PDE, choose the smooth physical state
    %
    %       u(s) = sin(pi*s),       v(s) = u_ss(s).
    %
    % Since u_t = u_ss + 5*u + u^2 and v=u_ss, differentiation gives
    %
    % v_t(s) = (pi^4-5*pi^2)*sin(pi*s) + 2*pi^2*cos(2*pi*s).
    %
    % Thus T*v_t=f(v), as required by the polynomial PIE, without having
    % to introduce a separate numerical inverse of T in this regression.
    v_test  = @(s)-pi^2*sin(pi*s);
    vt_test = @(s)(pi^4-5*pi^2)*sin(pi*s) + 2*pi^2*cos(2*pi*s);

    % Use identical quadrature for V3(v+h*vt), V3(v-h*vt), and dV3(v).
    % The centered difference is second-order in h.
    quad_order = 8;
    h_fd       = 5e-4;
    V3_plus = evaluate_fdp_numeric(V3,@(s)v_test(s)+h_fd*vt_test(s), ...
        dvar_names,dvar_values,quad_order);
    V3_minus = evaluate_fdp_numeric(V3,@(s)v_test(s)-h_fd*vt_test(s), ...
        dvar_names,dvar_values,quad_order);
    dV3_finite_difference = (V3_plus-V3_minus)/(2*h_fd);
    dV3_numeric = evaluate_fdp_numeric(dV3,v_test, ...
        dvar_names,dvar_values,quad_order);

    fd_scale = max([1,abs(dV3_numeric),abs(dV3_finite_difference)]);
    fd_error = abs(dV3_finite_difference-dV3_numeric)/fd_scale;
    assert(fd_error < 5e-4, ...
        'V3_Liediff_test:FiniteDifferenceMismatch', ...
        ['The d=%d analytical Lie derivative differs from the centered ', ...
         'finite difference by %.3e (relative scaled error).'],d,fd_error);

    fprintf(['     Passed d=%d, opdeg=%d: ', ...
             'deg(V3)={%s}, deg(dV3)={%s}, symmetry error=%.3e, ', ...
             'FD error=%.3e\n'], ...
        d,opdeg,num2str(V3_degs_actual'),num2str(dV3_degs_actual'), ...
        symmetry_error,fd_error);
end

fprintf([' --- All V3_Liediff construction, degree, analytical, ', ...
         'symmetry, and numerical regression tests passed ---\n\n']);


function [Zx,ZTx] = build_regression_factors(PIE,d,opdeg)
    % Independently construct U*x and U*T*x for the Lem. 13 regressions.

    pvar s s_dum
    Zmon     = monomials([s,s_dum],0:opdeg);
    Zop      = opvar();
    Zop.var1 = s;
    Zop.var2 = s_dum;
    Zop.I    = PIE.dom;

    if d == 1
        Zmon0    = monomials(s,0:opdeg);
        Zop.R.R0 = [Zmon0;0*Zmon;0*Zmon];
        Zop.R.R1 = [0*Zmon0;Zmon;0*Zmon];
        Zop.R.R2 = [0*Zmon0;0*Zmon;Zmon];
    else
        Zop.R.R0 = [0*Zmon;0*Zmon];
        Zop.R.R1 = [Zmon;0*Zmon];
        Zop.R.R2 = [0*Zmon;Zmon];
    end

    Z    = dopvar2ndopvar(Zop);
    Zx   = Z*PIE.f.vartab;
    ZTx  = (Z*PIE.T)*PIE.f.vartab;
end


function dF = substitute_pie_rhs(F,factor_idx,f)
    % Replace one differentiated state factor by sum_l C_l*x^(otimes l).

    dF = 0;
    for l = 1:size(f.degmat,1)
        f_l           = f;
        f_l.degmat    = f.degmat(l,:);
        f_l.C.ops     = f.C.ops(:,l);
        f_l.C.depmat2 = f.C.depmat2(l,:);
        dF = dF + subs(F,factor_idx,f_l);
    end
end


function verify_fdp_equality(actual,expected,d)
    % Canonicalize the two FDPs and verify equality of every kernel.

    residual = piesos_combine_terms(actual-expected);
    max_residual = 0;

    for term_idx = 1:numel(residual.C.ops)
        Kop = residual.C.ops{term_idx};
        if isempty(Kop)
            continue
        elseif isa(Kop,'intop')
            params = Kop.params;
        else
            params = Kop;
        end

        if isa(params,'dpvar')
            coeffs = params.C;
        elseif isa(params,'polynomial')
            coeffs = params.coefficient;
        elseif isnumeric(params)
            coeffs = params;
        else
            error('V3_Liediff_test:UnsupportedCoefficientClass', ...
                'Cannot compare FDP coefficients of class "%s".',class(params));
        end

        if ~isempty(coeffs)
            max_residual = max(max_residual,max(abs(full(coeffs(:)))));
        end
    end

    assert(max_residual < 1e-10, ...
        'V3_Liediff_test:IncorrectAnalyticalLieDerivative', ...
        ['The computed dV3 does not match the explicit d=%d ', ...
         'analytical expression from Lem. 13.'],d);
end


function [dvar_names,dvar_values] = symmetry_feasible_point(prog)
    % Construct a deterministic nonzero solution of the homogeneous
    % equality constraints representing Qhat_ij^*=Qhat_ji.

    dvar_names = prog.decvartable(:);
    ndvars = numel(dvar_names);
    [Aeq,beq] = collect_equality_system(prog);

    assert(norm(beq,inf)<1e-12, ...
        'V3_Liediff_test:NonhomogeneousSymmetryConstraint', ...
        'The Gram-block symmetry constraints were expected to be homogeneous.');

    % Project a reproducible trial vector onto null(Aeq). Using a dense
    % trigonometric sequence avoids dependence on MATLAB's random stream.
    trial = sin(sqrt(2)*(1:ndvars)') + cos(sqrt(3)*(1:ndvars)');
    if isempty(Aeq)
        dvar_values = trial;
    else
        dvar_values = trial-lsqminnorm(Aeq,Aeq*trial,1e-12);
    end

    max_value = max(abs(dvar_values));
    assert(max_value>1e-10, ...
        'V3_Liediff_test:TrivialSymmetryPoint', ...
        'The symmetry constraints produced only a numerically zero test point.');
    dvar_values = dvar_values/max_value;

    symmetry_residual = norm(Aeq*dvar_values-beq,inf);
    assert(symmetry_residual<1e-9, ...
        'V3_Liediff_test:InfeasibleSymmetryPoint', ...
        'The numerical Gram point violates symmetry by %.3e.',symmetry_residual);
end


function relative_residual = verify_symmetry_is_enforced(prog,x,residual,i,j)
    % Verify that RESIDUAL=0 is implied by the equality constraints in PROG.

    [Aeq,beq] = collect_equality_system(prog);

    % Convert the independently reconstructed operator identity into its
    % coefficient equations using a temporary program with the same ordered
    % decision-variable table.
    comparison_prog = piesos_program(x,dpvar(prog.decvartable));
    comparison_prog = piesos_eq(comparison_prog,residual);
    [Asym,bsym] = collect_equality_system(comparison_prog);

    assert(~isempty(Asym), ...
        'V3_Liediff_test:EmptySymmetryConstraint', ...
        'The reconstructed symmetry identity for block (%d,%d) is empty.',i,j);
    assert(norm(beq,inf)<1e-12 && norm(bsym,inf)<1e-12, ...
        'V3_Liediff_test:NonhomogeneousSymmetryConstraint', ...
        'The block (%d,%d) symmetry identity must be homogeneous.',i,j);

    % Rows of Asym must belong to row(Aeq). Equivalently, all solutions of
    % Aeq*q=0 must satisfy Asym*q=0. Solve all row-space projection problems
    % together and check the relative residual.
    row_coefficients = lsqminnorm(Aeq',Asym',1e-12);
    row_residual = Aeq'*row_coefficients-Asym';
    relative_residual = norm(row_residual,'fro')/max(1,norm(Asym,'fro'));
    assert(relative_residual<1e-9, ...
        'V3_Liediff_test:SymmetryNotEnforced', ...
        ['The equality constraints do not enforce ', ...
         'Qhat_{%d,%d}^*=Qhat_{%d,%d}; row-space residual %.3e.'], ...
        i,j,j,i,relative_residual);
end


function [Aeq,beq] = collect_equality_system(prog)
    % Stack the SOSTOOLS equality equations as Aeq*d=beq.

    ndvars = numel(prog.decvartable);
    Aeq = sparse(0,ndvars);
    beq = zeros(0,1);
    for constraint_idx = 1:prog.expr.num
        if strcmp(prog.expr.type{constraint_idx},'eq')
            Aeq = [Aeq;prog.expr.At{constraint_idx}']; %#ok<AGROW>
            beq = [beq;prog.expr.b{constraint_idx}]; %#ok<AGROW>
        end
    end
end


function value = evaluate_fdp_numeric(F,state_fcn,dvar_names,dvar_values,nquad)
    % Numerically evaluate a scalar, single-state FDP using Gauss quadrature.

    assert(isequal(F.matdim,[1,1]) && size(F.degmat,2)==1, ...
        'V3_Liediff_test:UnsupportedNumericalFDP', ...
        'The numerical regression supports scalar, single-state FDPs only.');

    value = 0;
    for mon_idx = 1:size(F.degmat,1)
        degree = F.degmat(mon_idx,1);
        Kop = F.C.ops{mon_idx};
        if isempty(Kop)
            continue
        end
        assert(isa(Kop,'intop') && isequal(Kop.matdim,[1,1]), ...
            'V3_Liediff_test:UnsupportedNumericalCoefficient', ...
            'Each nonconstant FDP coefficient must be a scalar intop.');

        params = evaluate_dpvar_decisions(Kop.params,dvar_names,dvar_values);
        for region_idx = 1:size(Kop.omat,1)
            kernel = params(1,region_idx);
            value = value + integrate_intop_region(kernel,Kop, ...
                Kop.omat(region_idx,:),degree,state_fcn,nquad);
        end
    end
end


function value = integrate_intop_region(kernel,Kop,order,degree,state_fcn,nquad)
    % Integrate one intop kernel over its ordered simplex. An all-zero
    % ordering denotes the quadratic multiplier term supported by intop.

    [nodes,weights] = gauss_legendre_rule(nquad);
    dom = Kop.dom(1,:);

    if all(order==0)
        assert(degree==2, ...
            'V3_Liediff_test:UnexpectedMultiplierDegree', ...
            'An intop multiplier term was encountered outside degree two.');
        spatial_nodes = (dom(2)-dom(1))*nodes/2 + sum(dom)/2;
        spatial_weights = (dom(2)-dom(1))*weights/2;
        points = repmat(spatial_nodes',degree,1);
        kernel_values = evaluate_kernel(kernel,Kop.pvarname,points);
        state_values = state_fcn(spatial_nodes').^degree;
        value = sum(spatial_weights'.*kernel_values.*state_values);
        return
    end

    assert(numel(order)==degree && ...
           isequal(sort(order),1:degree), ...
        'V3_Liediff_test:InvalidIntegrationOrder', ...
        'Each non-multiplier intop ordering must be a permutation.');

    % Tensor Gauss rule on [0,1]^degree.
    unit_nodes = (nodes+1)/2;
    unit_weights = weights/2;
    node_grids = cell(1,degree);
    weight_grids = cell(1,degree);
    [node_grids{:}] = ndgrid(unit_nodes);
    [weight_grids{:}] = ndgrid(unit_weights);
    npoints = numel(node_grids{1});
    cube_nodes = zeros(degree,npoints);
    cube_weights = ones(1,npoints);
    for idx = 1:degree
        cube_nodes(idx,:) = node_grids{idx}(:)';
        cube_weights = cube_weights.*weight_grids{idx}(:)';
    end

    % Map the cube recursively to a<=y1<=...<=y_degree<=b.
    ordered_nodes = zeros(degree,npoints);
    jacobian = ones(1,npoints);
    lower_limit = dom(1)*ones(1,npoints);
    for idx = 1:degree
        interval_length = dom(2)-lower_limit;
        ordered_nodes(idx,:) = lower_limit + ...
            interval_length.*cube_nodes(idx,:);
        jacobian = jacobian.*interval_length;
        lower_limit = ordered_nodes(idx,:);
    end

    % order=[k1,...,kd] represents t_k1<=...<=t_kd.
    points = zeros(degree,npoints);
    for idx = 1:degree
        points(order(idx),:) = ordered_nodes(idx,:);
    end

    kernel_values = evaluate_kernel(kernel,Kop.pvarname,points);
    state_values = prod(state_fcn(points),1);
    value = sum(cube_weights.*jacobian.*kernel_values.*state_values);
end


function values = evaluate_kernel(kernel,pvar_names,points)
    % Evaluate a scalar polynomial kernel at several spatial points.

    if isnumeric(kernel)
        values = kernel*ones(1,size(points,2));
        return
    end

    assert(isa(kernel,'polynomial'), ...
        'V3_Liediff_test:NonPolynomialKernel', ...
        'Numerical FDP evaluation requires polynomial kernels.');
    if isempty(kernel.varname)
        values = double(kernel)*ones(1,size(points,2));
        return
    end

    [is_present,var_indices] = ismember(kernel.varname,pvar_names);
    assert(all(is_present), ...
        'V3_Liediff_test:UnknownKernelVariable', ...
        'A kernel variable is absent from the intop spatial variables.');
    values = double(subs(kernel,kernel.varname,points(var_indices,:)));
    values = reshape(values,1,[]);
end


function P = evaluate_dpvar_decisions(D,dvar_names,dvar_values)
    % Substitute a numerical decision vector into a dpvar matrix while
    % retaining its polynomial dependence on the spatial variables.

    if ~isa(D,'dpvar')
        P = D;
        return
    end

    m = D.matdim(1);
    n = D.matdim(2);
    nmon = size(D.degmat,1);
    [is_present,dvar_indices] = ismember(D.dvarname,dvar_names);
    assert(all(is_present), ...
        'V3_Liediff_test:MissingDecisionVariable', ...
        'An FDP decision variable is absent from the test decision vector.');

    local_values = dvar_values(dvar_indices);
    Cnumeric = kron(speye(m),sparse([1;local_values(:)]'))*D.C;
    if nmon==0 || isempty(D.varname)
        P = reshape(full(Cnumeric),m,n);
        return
    end

    [row_idx,col_idx,coeff_values] = find(Cnumeric);
    matrix_col = floor((col_idx-1)/nmon);
    monomial_idx = col_idx-matrix_col*nmon;
    coefficients = sparse(monomial_idx,row_idx+m*matrix_col, ...
        coeff_values,nmon,m*n);
    P = polynomial(coefficients,D.degmat,D.varname,[m,n]);
end


function [nodes,weights] = gauss_legendre_rule(order)
    % Golub-Welsch construction of an order-point rule on [-1,1].

    indices = (1:order-1)';
    off_diagonal = indices./sqrt(4*indices.^2-1);
    Jacobi = diag(off_diagonal,1)+diag(off_diagonal,-1);
    [vectors,values] = eig(Jacobi);
    [nodes,permutation] = sort(diag(values));
    weights = 2*vectors(1,permutation)'.^2;
end

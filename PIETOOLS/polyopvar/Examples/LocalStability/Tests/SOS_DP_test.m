%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - SOS_DP_test.m
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

% Test the SOS distributed polynomials returned by SOS_DP.m against Def. 9
% of DSOS.pdf. The test verifies the U-hat basis, tensor-power dimensions,
% PSD Gram declaration, complete quadratic form, and nonnegativity of
% z_d(v)'*Q*z_d(v) for representative states and nonzero PSD matrices Q.

clearvars; clear stateNameGenerator; close all; clc;

% Add PIETOOLS to the path when this test is run directly from its Tests
% directory or through MATLAB's run command.
test_dir      = fileparts(mfilename('fullpath'));
local_dir     = fileparts(test_dir);
pietools_root = fileparts(fileparts(fileparts(local_dir)));
addpath(genpath(pietools_root));


%% 1. Construct a scalar distributed state.

% Use the scalar fundamental state of the Fisher PIE, matching the state
% space assumed in Def. 9 and in the LocalStability examples.
pvar s t
dom = [0,1];
u   = pde_var(s,dom);
PDE = [diff(u,t)==diff(u,s,2) + 5*u + u^2;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];

PIE = convert(PDE);
x   = PIE.f.vartab;


%% 2. Test low-order distributed and kernel monomial degrees.

% The d=2 case exercises the recursive distributed-monomial basis. The two
% d=1 cases check both constant and nonconstant U_bar_d(s,theta) bases.
test_cases = [1,0;
              1,1;
              2,0;
              2,1];

state_fcns = {@(theta)sin(pi*theta), ...
              @(theta)theta.*(1-theta), ...
              @(theta)cos(2*pi*theta)+0.25*sin(pi*theta)};

fprintf('\n --- Running SOS_DP Definition 9 and nonnegativity tests ---\n');

for test_idx = 1:size(test_cases,1)
    d     = test_cases(test_idx,1);
    opdeg = test_cases(test_idx,2);

    clear stateNameGenerator
    prog = piesos_program(x);
    [prog,DP,Pcell,Zs] = SOS_DP(prog,d,opdeg,x,dom);

    %% Definition 9 basis and block dimensions.

    mu = (opdeg+1)*(opdeg+2)/2;
    k  = 2*mu;
    block_dims = k.^(1:d);

    assert(isa(DP,'polyopvar') && isequal(size(DP),[1,1]), ...
        'SOS_DP_test:InvalidPolynomial', ...
        'SOS_DP must return a scalar polyopvar object.');
    assert(isequal(size(Pcell),[d,d]) && isequal(size(Zs),[d,1]), ...
        'SOS_DP_test:InvalidBlockStructure', ...
        'Pcell and Zs do not have the block structure required by Def. 9.');

    for i = 1:d
        assert(isequal(Zs{i}.degmat,i) && size(Zs{i},1)==block_dims(i), ...
            'SOS_DP_test:IncorrectBasisPower', ...
            ['Zs{%d} must represent U-hat^(otimes %d)x^(otimes %d) ', ...
             'with %d components.'],i,i,i,block_dims(i));
        for j = 1:d
            assert(isequal(size(Pcell{i,j}),[block_dims(i),block_dims(j)]), ...
                'SOS_DP_test:IncorrectGramBlockSize', ...
                'Gram block (%d,%d) has an incorrect Def. 9 dimension.',i,j);
        end
    end

    expected_degrees = (2:2*d)';
    actual_degrees = unique(sum(DP.degmat,2));
    assert(isequal(actual_degrees(:),expected_degrees), ...
        'SOS_DP_test:IncorrectDistributedDegrees', ...
        'SOS_DP does not contain the distributed degrees required by Z_d.');

    %% PSD Gram declaration and complete Def. 9 quadratic form.

    P = concatenate_blocks(Pcell);
    gram_dim = sum(block_dims);
    assert(isa(P,'dpvar') && isequal(size(P),[gram_dim,gram_dim]), ...
        'SOS_DP_test:IncorrectGramDimension', ...
        'The global Gram matrix has an incorrect dimension.');
    verify_dpvar_zero(P-P', ...
        'SOS_DP_test:NonsymmetricGram', ...
        'The Gram matrix declared by SOS_DP is not symmetric.');
    assert(strcmp(prog.var.type{end},'sos'), ...
        'SOS_DP_test:GramNotPSD', ...
        'The Gram matrix must be declared as a PSD SOSTOOLS variable.');

    % Independently expand the full i,j block sum in (15)-(16). SOS_DP uses
    % symmetry to evaluate only one triangle; both forms must be identical.
    DP_expected = 0;
    for i = 1:d
        for j = 1:d
            DP_expected = DP_expected + ...
                innerprod_v2(Zs{i},Zs{j},Pcell{i,j});
        end
    end
    verify_fdp_equality(DP,DP_expected, ...
        'SOS_DP_test:IncorrectQuadraticForm', ...
        'SOS_DP does not equal the complete Def. 9 block quadratic form.');

    %% Numerical nonnegativity regression.

    % Every PSD Q admits Q=A'*A, and therefore the Def. 9 integrand is
    % z_d(v)'*Q*z_d(v)=||A*z_d(v)||^2>=0. Test the identity matrix, a
    % rank-one matrix, and a deterministic dense positive-definite matrix.
    direction = sin(sqrt(2)*(1:gram_dim)');
    dense_factor = reshape(sin(sqrt(3)*(1:gram_dim^2)),gram_dim,gram_dim);
    Qtests = {eye(gram_dim), ...
              direction*direction', ...
              dense_factor'*dense_factor+0.1*eye(gram_dim)};

    minimum_value = inf;
    for Q_idx = 1:numel(Qtests)
        Q = Qtests{Q_idx};
        assert(min(eig((Q+Q')/2))>-1e-10, ...
            'SOS_DP_test:InvalidNumericalGram', ...
            'A numerical regression Gram matrix is not PSD.');

        for state_idx = 1:numel(state_fcns)
            q_value = evaluate_def9_quadratic( ...
                state_fcns{state_idx},Q,d,opdeg,dom,10);
            minimum_value = min(minimum_value,q_value);
            assert(q_value>=-1e-9, ...
                'SOS_DP_test:NegativeSOSPolynomial', ...
                ['The Def. 9 quadratic form is negative for d=%d, ', ...
                 'opdeg=%d: q(v)=%.3e.'],d,opdeg,q_value);
        end
    end

    fprintf(['     Passed d=%d, opdeg=%d: Gram dimension=%d, ', ...
             'minimum sampled q(v)=%.3e\n'], ...
        d,opdeg,gram_dim,minimum_value);
end

fprintf([' --- All SOS_DP Definition 9 and nonnegativity tests ', ...
         'passed ---\n\n']);


function P = concatenate_blocks(Pcell)
    % Reassemble the global Gram matrix from its Pcell blocks.

    block_rows = cell(size(Pcell,1),1);
    for i = 1:size(Pcell,1)
        block_rows{i} = horzcat(Pcell{i,:});
    end
    P = vertcat(block_rows{:});
end


function verify_dpvar_zero(P,error_id,error_message)
    % Verify that every affine coefficient in a dpvar matrix is zero.

    assert(isa(P,'dpvar') && ...
           (isempty(P.C) || max(abs(full(P.C(:))))<1e-12), ...
        error_id,error_message);
end


function verify_fdp_equality(actual,expected,error_id,error_message)
    % Canonicalize and compare two decision-variable-valued FDPs.

    residual = piesos_combine_terms(actual-expected);
    max_residual = 0;
    for term_idx = 1:numel(residual.C.ops)
        coefficient = residual.C.ops{term_idx};
        if isempty(coefficient)
            continue
        elseif isa(coefficient,'intop')
            coefficient = coefficient.params;
        end

        if isa(coefficient,'dpvar')
            values = coefficient.C;
        elseif isa(coefficient,'polynomial')
            values = coefficient.coefficient;
        elseif isnumeric(coefficient)
            values = coefficient;
        else
            error('SOS_DP_test:UnsupportedCoefficient', ...
                'Cannot compare an FDP coefficient of class "%s".', ...
                class(coefficient));
        end

        if ~isempty(values)
            max_residual = max(max_residual,max(abs(full(values(:)))));
        end
    end
    assert(max_residual<1e-10,error_id,error_message);
end


function q_value = evaluate_def9_quadratic(state_fcn,Q,d,opdeg,dom,nquad)
    % Evaluate integral z_d(v)(s)'*Q*z_d(v)(s) ds from Def. 9 directly.

    [outer_nodes,outer_weights] = gauss_legendre_rule(nquad,dom(1),dom(2));
    q_values = zeros(size(outer_nodes));
    for node_idx = 1:numel(outer_nodes)
        Uv = evaluate_Uhat(state_fcn,outer_nodes(node_idx),opdeg,dom,nquad);
        Zd = Uv;
        Upower = Uv;
        for degree_idx = 2:d
            Upower = kron(Upower,Uv);
            Zd = [Zd;Upower]; %#ok<AGROW>
        end
        q_values(node_idx) = Zd'*Q*Zd;
    end
    q_value = sum(outer_weights.*q_values);
end


function Uv = evaluate_Uhat(state_fcn,sval,opdeg,dom,nquad)
    % Evaluate the two integral components of U-hat in Def. 9 at s=sval.

    pvar s theta
    Umon = monomials([s,theta],0:opdeg);

    [lower_nodes,lower_weights] = gauss_legendre_rule( ...
        nquad,dom(1),sval);
    [upper_nodes,upper_weights] = gauss_legendre_rule( ...
        nquad,sval,dom(2));

    lower_basis = double(subs(Umon,[s;theta], ...
        [sval*ones(1,numel(lower_nodes));lower_nodes']));
    upper_basis = double(subs(Umon,[s;theta], ...
        [sval*ones(1,numel(upper_nodes));upper_nodes']));

    lower_value = lower_basis*(lower_weights.*state_fcn(lower_nodes));
    upper_value = upper_basis*(upper_weights.*state_fcn(upper_nodes));
    Uv = [lower_value;upper_value];
end


function [nodes,weights] = gauss_legendre_rule(order,a,b)
    % Golub-Welsch construction of an order-point rule on [a,b].

    indices = (1:order-1)';
    off_diagonal = indices./sqrt(4*indices.^2-1);
    Jacobi = diag(off_diagonal,1)+diag(off_diagonal,-1);
    [vectors,values] = eig(Jacobi);
    [unit_nodes,permutation] = sort(diag(values));
    unit_weights = 2*vectors(1,permutation)'.^2;
    nodes = (b-a)*unit_nodes/2+(a+b)/2;
    weights = (b-a)*unit_weights/2;
end

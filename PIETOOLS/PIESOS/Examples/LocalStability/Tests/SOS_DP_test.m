%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - SOS_DP_test
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
% DJ, 10/15/2026: Updated to use new polyopvar functions;

% Test the SOS distributed polynomials returned by SOS_DP.m against Def. 9
% of DSOS.pdf. The test verifies the block dimensions of the basis and
% Gram matrix, the PSD Gram declaration, the complete quadratic form, and
% the SOS property by feasibility programs.

clear stateNameGenerator;
all_pass = true;


%% 1. Construct a scalar distributed state.

% Use the scalar fundamental state of the Fisher PIE, matching the state
% space assumed in Def. 9 and in the LocalStability examples.
pvar s t
dom = [0,1];
u   = pde_var(s,dom);
PDE = [diff(u,t)==diff(u,s,2) + 5*u + u^2;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];

PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
x   = PIE.x;
Top = PIE.T;
Top_o = sopvar2opvar(Top.ops{1}.ops{1,1});
Rop = tensopmat(sopvar2tensopvar(opvar2sopvar(diff(Top_o,Top_o.var1,1,'pure'))));
uu  = Top*x;    uus = Rop*x;


%% 2. Test low-order distributed and kernel monomial degrees.

% The d=2 case exercises the recursive distributed-monomial basis. The two
% d=1 cases check both constant and nonconstant U_bar_d(s,theta) bases.
test_cases = [1,0;
              1,1;
              2,0;
              2,1];

fprintf('\n --- Running SOS_DP Definition 9 tests ---\n');

for test_idx = 1:size(test_cases,1)
    d     = test_cases(test_idx,1);
    opdeg = test_cases(test_idx,2);

    clear stateNameGenerator
    prog = piesos_program(x);
    [prog,DP,Pcell,Zs] = SOS_DP(prog,d,opdeg,x,dom);

    %% Definition 9 basis and block dimensions.

    % Components of U^i*x^i in SOS_DP: for i=1, 2 integral terms with the monomials in (s,theta)
    % up to degree opdeg in each (m1+1)^2 monomials, m1=opdeg; for i>1 the first factor has the monomials
    % up to degree floor(opdeg/i) in each (m+1)^2, and the other i-1 factors 2*(m+1) components each
    block_dims = zeros(1,d);
    for i = 1:d
        m = floor(opdeg/i);
        if i==1
            block_dims(i) = 2*(opdeg+1)^2;
        else
            block_dims(i) = 2*(m+1)^2*(2*(m+1))^(i-1);
        end
    end

    assert(isa(DP,'polyopvar') && isequal(size(DP),[1,1]), ...
        'SOS_DP_test:InvalidPolynomial', ...
        'SOS_DP must return a scalar polyopvar object.');
    assert(isequal(size(Pcell),[d,d]) && isequal(size(Zs),[d,1]), ...
        'SOS_DP_test:InvalidBlockStructure', ...
        'Pcell and Zs do not have the block structure required by Def. 9.');

    for i = 1:d
        assert(isequal(full(Zs{i}.degmat),i) && size(Zs{i},1)==block_dims(i), ...
            'SOS_DP_test:IncorrectBasisPower', ...
            ['Zs{%d} must represent U-hat^(otimes %d)x^(otimes %d) ', ...
             'with %d components, has %d.'],i,i,i,block_dims(i),size(Zs{i},1));
        for j = 1:d
            assert(isequal(size(Pcell{i,j}),[block_dims(i),block_dims(j)]), ...
                'SOS_DP_test:IncorrectGramBlockSize', ...
                'Gram block (%d,%d) has an incorrect Def. 9 dimension.',i,j);
        end
    end

    expected_degrees = (2:2*d)';
    actual_degrees = unique(sum(DP.degmat,2));
    assert(isequal(full(actual_degrees(:)),expected_degrees), ...
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

    fprintf('     Passed d=%d, opdeg=%d: Gram dimension=%d\n',d,opdeg,gram_dim);
end


%% 3. SOS property: feasibility programs DP = target.

% Targets that are SOS by construction (feasible) or not (infeasible).
% SOS_DP has no multiplier term in its basis operators, so targets without
% delta terms (the kernel terms int (Tx)^2 and int (Rx)^2) are used.
target_k = innerprod(uu,uu) + innerprod(uus,uus);
w        = times(uu,uus);
target_4 = target_k + 3*innerprod(w,w);
% weighted: int s(1-s)*(Tx)^2 ds, of the form g*P for psatz=0 (the weight g=(s-a)*(b-s) vanishes
% at the boundary, so that psatz=0 cannot represent the unweighted target_k)
target_g = innerprod(uu,s*(1-s),uu);

% each row: d, opdeg, opts, target, expected feasibility
cases = {
 'kernel terms',    1, 1, [],                    target_k,   true;
 '-kernel terms',   1, 1, [],                    -target_k,  false;
 'degree 4, d=2',   2, 2, [],                    target_4,   true;
 'degree 4, d=1',   1, 1, [],                    target_4,   false;
 'psatz=0 (g*P)',   1, 1, struct('psatz',0),    target_g,   true;
 'psatz=0, no weight', 1, 1, struct('psatz',0),   target_k,   false;
 'psatz=1 (P0+g*P1)',1,1, struct('psatz',1),    target_k,   true;
};
nfail = 0;
for k=1:size(cases,1)
    clear stateNameGenerator
    prog = piesos_program(x);
    if isempty(cases{k,4})
        [prog,DPk] = SOS_DP(prog,cases{k,2},cases{k,3},x,dom);
    else
        % opts must have the field psatz, if specified
        [prog,DPk] = SOS_DP(prog,cases{k,2},cases{k,3},x,dom,cases{k,4});
    end
    prog = piesos_eq(prog,DPk-cases{k,5});
    sol_opts = struct('solver','mosek','simplify',true);
    evalc('prog_sol = piesos_solve(prog,sol_opts);');
    info = prog_sol.solinfo.info;
    feasible = ~info.pinf && ~info.dinf && ~info.numerr;
    status = 'ok';
    if feasible~=cases{k,6}
        status = 'FAIL';    nfail = nfail+1;
    end
    fprintf('     %-4s %-20s expected feasible=%d: pinf=%d dinf=%d numerr=%d\n', ...
        status,cases{k,1},cases{k,6},info.pinf,info.dinf,info.numerr);
end
assert(nfail==0,'SOS_DP_test:Infeasible','%d of the SOS_DP feasibility cases failed.',nfail);

fprintf([' --- All SOS_DP Definition 9 and SOS tests ', ...
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

    [is_equal,~] = polyopvar_equal(actual,expected,1e-10);
    assert(is_equal,error_id,error_message);
end
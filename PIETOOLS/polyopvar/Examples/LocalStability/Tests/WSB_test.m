%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - WSB_test.m
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

% Test the weighted Sobolev bound and local-domain polynomial used in
% LocalStability.m against Theorem 1 of Auto.pdf. In the notation of the
% theorem,
%
%   B(v)   = sum_{j=0}^n alpha_j vec(R_j^(otimes 2)) v^(otimes 2),
%   g_r(v) = r^2-B(v),                 R_j = partial_s^j o T.
%
% LocalStability.m stores Weighted_Sobolev_Ball(0,...) in the variable
% named "bound". That object equals -B(v), so this test checks both the
% positive theorem quantity B(v) and the negated convention used by the
% lower- and upper-bound inequalities in LocalStability.m.

clearvars; clear stateNameGenerator; close all; clc;

% Add PIETOOLS to the path when this test is run directly from its Tests
% directory or by using MATLAB's run command.
test_dir      = fileparts(mfilename('fullpath'));
local_dir     = fileparts(test_dir);
pietools_root = fileparts(fileparts(fileparts(local_dir)));
addpath(genpath(pietools_root));


%% 1. Construct a second-order scalar PDE and its PIE representation.

% The Fisher equation has spatial order n=2, allowing the test to exercise
% R_0=T, R_1=partial_s o T, and R_2=partial_s^2 o T independently.
pvar s t
dom = [0,1];
u   = pde_var(s,dom);
PDE = [diff(u,t)==diff(u,s,2) + 5*u + u^2;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];

PIE = convert(PDE);
Top = PIE.T;
x   = PIE.f.vartab;


%% 2. Independently construct the operators R_j from Theorem 1.

n = 2;
Rops = cell(n+1,1);
Rops{1} = Top;                         % R_0 = T.

Top_opvar = ndopvar2dopvar(Top);
for j = 1:n
    % MATLAB cell index j+1 corresponds to the paper index j.
    Rops{j+1} = dopvar2ndopvar( ...
        diff(Top_opvar,Top_opvar.var1,j,'pure'));
end

% Form each unweighted theorem term
%
%       B_j(v) = vec(R_j^(otimes 2))v^(otimes 2)
%              = <R_j*v,R_j*v>_L2.
%
% Constructing these terms before calling Weighted_Sobolev_Ball keeps the
% expected expression independent of that function's loop and signs.
B_terms = cell(n+1,1);
for j = 0:n
    Rjx = Rops{j+1}*x;
    B_terms{j+1} = innerprod_v2(Rjx,Rjx);
end


%% 3. Check individual and mixed Sobolev weights.

r = 1.7;

% Unit weights isolate the alpha_0, alpha_1, and alpha_2 indexing. The
% mixed case checks accumulation of all terms and non-unit coefficients.
alpha_cases = [1,0,0;
               0,1,0;
               0,0,1;
               2,3,4];

fprintf('\n --- Running LocalStability weighted Sobolev tests ---\n');

for case_idx = 1:size(alpha_cases,1)
    alpha = alpha_cases(case_idx,:);

    % Positive bound B(v) appearing in Theorem 1.
    B_expected = 0;
    for j = 0:n
        B_expected = B_expected + alpha(j+1)*B_terms{j+1};
    end
    g_expected = r^2-B_expected;

    % These are exactly the two calls made in LocalStability.m.
    g_actual = Weighted_Sobolev_Ball(r,alpha,Top,x);
    bound_local = Weighted_Sobolev_Ball(0.0,alpha,Top,x);

    % The first identity checks g_r(v)=r^2-B(v). The second confirms that
    % the variable called "bound" in LocalStability is intentionally -B(v).
    verify_fdp_equality(g_actual,g_expected, ...
        'LocalStability_bound_test:IncorrectBall', ...
        'Weighted_Sobolev_Ball does not equal r^2-B(v).');
    verify_fdp_equality(-bound_local,B_expected, ...
        'LocalStability_bound_test:IncorrectBound', ...
        ['The variable bound stored in LocalStability does not equal ', ...
         '-sum_j alpha_j vec(R_j^(otimes 2))v^(otimes 2).']);

    % This equivalent identity checks the sign convention used directly in
    % LocalStability: g = r^2+bound_local because bound_local=-B.
    verify_fdp_equality(g_actual,r^2+bound_local, ...
        'LocalStability_bound_test:InconsistentSignConvention', ...
        'The LocalStability relation g=r^2+bound is not satisfied.');

    assert(isequal(unique(sum(B_expected.degmat,2)),2), ...
        'LocalStability_bound_test:IncorrectBoundDegree', ...
        'The theorem bound B(v) must have distributed degree two.');
    assert(isequal(unique(sum(g_actual.degmat,2)),[0;2]), ...
        'LocalStability_bound_test:IncorrectBallDegrees', ...
        'The local-domain polynomial g_r(v) must have degrees zero and two.');

    fprintf('     Passed alpha=[%s]\n',num2str(alpha));
end

fprintf([' --- All LocalStability weighted Sobolev bound and ball ', ...
         'tests passed ---\n\n']);


function verify_fdp_equality(actual,expected,error_id,error_message)
    % Verify exact equality of two FDPs after common-basis addition.

    residual = actual-expected;
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
            error('LocalStability_bound_test:UnsupportedCoefficient', ...
                'Cannot compare an FDP coefficient of class "%s".', ...
                class(coefficient));
        end

        if ~isempty(values)
            max_residual = max(max_residual,max(abs(full(values(:)))));
        end
    end

    assert(max_residual<1e-12,error_id,error_message);
end

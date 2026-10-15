%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - WSB_test
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

% Test the weighted Sobolev bound and local-domain polynomial used in
% LocalStability.m against Theorem 1 of Auto.pdf. In the notation of the
% theorem,
%
%   B(v)   = sum_{j=0}^n alpha_j vec(R_j^(otimes 2)) v^(otimes 2),
%   g_r(v) = r^2-B(v),                 R_j = partial_s^j o T.
%
% LocalStability.m stores Weighted_Sobolev_Norm(...) in the variable
% named "bound", which equals B(v),
% and defines g = r^2 - bound.

clearvars; clear stateNameGenerator; close all; clc;


%% 1. Construct a second-order scalar PDE and its PIE representation.

% The Fisher equation has spatial order n=2, allowing the test to exercise
% R_0=T, R_1=partial_s o T, and R_2=partial_s^2 o T independently.
pvar s t
dom = [0,1];
u   = pde_var(s,dom);
PDE = [diff(u,t)==diff(u,s,2) + 5*u + u^2;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];

PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;
x   = PIE.x;


%% 2. Independently construct the operators R_j from Theorem 1.

n = 2;
Rops = cell(n+1,1);
Rops{1} = Top;                         % R_0 = T.

Top_opvar = sopvar2opvar(Top.ops{1}.ops{1,1});
for j = 1:n
    % MATLAB cell index j+1 corresponds to the paper index j.
    Rops{j+1} = tensopmat(sopvar2tensopvar(opvar2sopvar( ...
        diff(Top_opvar,Top_opvar.var1,j,'pure'))));
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
    B_terms{j+1} = innerprod(Rjx,Rjx);
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

    % These are exactly the calls made in LocalStability.m.
    bound_local = Weighted_Sobolev_Norm(alpha,Top,x);
    g_actual = r^2 - bound_local;

    % The first identity checks g_r(v)=r^2-B(v). The second confirms that
    % the variable called "bound" in LocalStability is B(v).
    verify_fdp_equality(g_actual,g_expected, ...
        'LocalStability_bound_test:IncorrectBall', ...
        'r^2 - Weighted_Sobolev_Norm does not equal r^2-B(v).');
    verify_fdp_equality(bound_local,B_expected, ...
        'LocalStability_bound_test:IncorrectBound', ...
        ['The variable bound stored in LocalStability does not equal ', ...
         'sum_j alpha_j vec(R_j^(otimes 2))v^(otimes 2).']);
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

    % The residual is canonicalized and its coefficients (kernels) read by
    % polyopvar_equal, which handles the 'intop' kernels of the new classes
    [is_equal,~] = polyopvar_equal(actual,expected,1e-12);
    assert(is_equal,error_id,error_message);
end

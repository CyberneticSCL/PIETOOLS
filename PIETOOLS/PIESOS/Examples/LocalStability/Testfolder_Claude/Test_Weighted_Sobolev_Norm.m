function Test_Weighted_Sobolev_Norm()
% TEST_WEIGHTED_SOBOLEV_NORM verifies Weighted_Sobolev_Norm(alpha,Top,x),
% the weighted Sobolev norm B(x) = sum_j alpha_j*||d^j/ds^j (T*x)||^2_{L2}
% of the PDE state u=T*x, for the PIE state x=u_ss, against an independent
% ground truth.
%
% Ground truth: for a PDE state u(s) satisfying the boundary conditions,
% the PIE state is x=u_ss, and T*x=u. B(x) is evaluated at x with
% eval_polyopvar (apply_op of the operators, not `subs`), and compared
% with the integrals of the squared derivatives of u, computed from the
% polynomial coefficients of u with polyder/polyint (MATLAB, independent
% of the PIETOOLS classes and of the derivative of the operator T that is
% under test, 'diff' for 'tensopmat'). The cases cover each of the derivatives j=0,1,2
% separately (the j=2 case has the multiplier term of d^2/ds^2 T, which
% equals the identity) and a mixed weight, for Dirichlet-Dirichlet and
% Neumann-Dirichlet boundary conditions (different operators T), and the
% errors for a map T that is not a single 'sopvar' operator.
%
% No SDPs are solved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_Weighted_Sobolev_Norm
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/20/2026: Weighted_Sobolev_Norm now differentiates T with
%                   'diff' for 'tensopmat': updated the description;

alpha_cases = [1,0,0; 0,1,0; 0,0,1; 2,3,4; 0,0,0];

% BC cases: name, BCs, polynomial coefficients (in s, highest power first) of
% test states u satisfying the BCs
% Dirichlet-Dirichlet: u(0)=u(1)=0; Neumann-Dirichlet: u_s(0)=0, u(1)=0
pvar s t
dom = [0,1];
bc_cases = struct('name',{'dirichlet','neumann_dirichlet'}, ...
    'u',{{[-1,1,0], [1,-2,1,0,0], [-1,3,-3,1,0,0,0]}, ...       % u = s(1-s), s^2(1-s)^2, s^3(1-s)^3
         {[-1,0,1], [1,0,-2,0,1], [-1,0,3,0,-3,0,1]}});         % u = (1-s^2)^k, k=1,2,3
nfail = 0;
for ib=1:numel(bc_cases)
    clear stateNameGenerator
    x = pde_var(s,dom);
    if ib==1
        PDE = [diff(x,t)==diff(x,s,2); subs(x,s,dom(1))==0; subs(x,s,dom(2))==0];
    else
        PDE = [diff(x,t)==diff(x,s,2); subs(diff(x,s),s,dom(1))==0; subs(x,s,dom(2))==0];
    end
    PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
    Top = PIE.T;
    xs = PIE.x;
    for ic=1:numel(bc_cases(ib).u)
        pu = bc_cases(ib).u{ic};                    % coefficients of u(s)
        % PIE state x = u_ss as a polynomial in s
        xval = polynomial(poly2pol(polyder(polyder(pu)),s));
        for ia=1:size(alpha_cases,1)
            alpha = alpha_cases(ia,:);
            B = Weighted_Sobolev_Norm(alpha,Top,xs);
            % Expected value: sum_j alpha_j*int_0^1 (u^{(j)})^2 ds
            Bexp = 0;
            for j=0:2
                dj = pu;
                for k=1:j,  dj = polyder(dj);   end
                sq = polyint(conv(dj,dj));
                Bexp = Bexp + alpha(j+1)*(polyval(sq,1)-polyval(sq,0));
            end
            if all(alpha==0)
                % B = 0 (the double 0, no terms)
                assert(isnumeric(B) && B==0,'Expected 0 for alpha=0.');
                Bval = 0;
            else
                assert(isa(B,'polyopvar') && isequal(size(B),[1,1]), ...
                    'Weighted_Sobolev_Norm must return a scalar polyopvar.');
                assert(all(sum(B.degmat,2)==2), ...
                    'The weighted Sobolev norm must be homogeneous of degree 2.');
                Bval = double(eval_V(B,xval));
            end
            err = abs(Bval-Bexp);
            if err>1e-10*max(1,abs(Bexp))
                nfail = nfail+1;
                fprintf('FAIL %s u#%d alpha=[%s]: B=%.12g, expected %.12g\n', ...
                    bc_cases(ib).name,ic,num2str(alpha),Bval,Bexp);
            end
        end
    end
    fprintf('Passed BC case: %s\n',bc_cases(ib).name);
end
assert(nfail==0,'%d weighted Sobolev norm cases failed.',nfail);

% Errors for inappropriate maps
threw = false;
try
    Weighted_Sobolev_Norm([1,0,0],1,xs);
catch
    threw = true;
end
assert(threw,'An error was expected for a map T that is not a tensopmat.');

disp('Test_Weighted_Sobolev_Norm: all tests passed.')
end


function p = poly2pol(c,s)
% P = POLY2POL(C,S) polynomial in the variable S with coefficients C of
% decreasing power, as for polyval.
n = numel(c);
p = 0*s;
for k=1:n
    p = p + c(k)*s^(n-k);
end
end

function Y = eval_V(V,xval)
% Y = EVAL_V(V,XVAL) evaluates the scalar functional V (polyopvar of 'intop'
% coefficients, in a single scalar state depending on the variable s) at the 
% state x=XVAL, polynomial in s. For a monomial of d factors, the dummy
% variables of the factors are s_1,...,s_d, so that the input to the
% operator is kron(x(s_1),...,x(s_d)) (as in apply_op of 'intop' and the 
% evaluators of Test_piesos_LF).
Zf = full(V.Z);
Y = 0;
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    d = nnz(Zf(k,:));
    v = 1;
    for j=1:d
        xj = xval;
        if d>1,     xj = subs(xval,pvar('s'),pvar(sprintf('s_%d',j)));     end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end
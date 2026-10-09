function [Cmono,relrms,d] = chebfit_mono(x,Y,ab,deg0,degmax,tol,x2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CMONO,RELRMS,D] = CHEBFIT_MONO(X,Y,AB,DEG0,DEGMAX,TOL[,X2]) least-squares
% polynomial fit of sampled functions in the Chebyshev basis on AB = [a,b],
% the degree raised from DEG0 in steps of 2 until the relative RMS residual
% is below TOL or D = DEGMAX, returned as MONOMIAL coefficients.
%   one variable:   Y is numel(X) x nc, Y(i,:) = f(X(i));
%                   CMONO is (d+1) x nc, f(s) ~ sum_r CMONO(r+1,:) s^r;
%   two variables:  X2 given, Y is numel(X)*numel(X2) x nc with row index
%                   i + (j-1)*numel(X) for (X(i),X2(j)); CMONO is (d+1)^2 x
%                   nc with row index (r-1)(d+1)+u for s^(r-1) t^(u-1).
% The Chebyshev basis keeps the least squares well conditioned at degree
% 16 on 101 points, where a monomial Vandermonde does not; the conversion
% to monomials is exact arithmetic on the recurrence.
%
% See also SOPVAR/INV, GETCONTROLLER_DIRECT_SOP, GK_GRID.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - chebfit_mono
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
% Initial coding MMP, 10/09/2026: the fit part of @sopvar/inv (same date),
%                moved here so that the gain construction can share it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = ab(1);  b = ab(2);
two = nargin>=7 && ~isempty(x2);
xs = (2*(x(:) - a)/(b - a)) - 1;
if two,     xt = (2*(x2(:) - a)/(b - a)) - 1;   end
d = deg0;
while true
    T1 = chebV(xs,d);
    if two
        T2 = chebV(xt,d);
        Phi = design2(T1,T2,d);
    else
        Phi = T1;
    end
    c = Phi\Y;
    res = Phi*c - Y;
    relrms = sqrt(mean(res(:).^2))/max(1e-300,sqrt(mean(Y(:).^2)));
    if ~any(Y(:)),  relrms = 0;     end
    if relrms<=tol || d>=degmax,    break,  end
    d = min(d+2,degmax);
end
Mc = cheb2mono(d,2/(b-a),-(a+b)/(b-a));
if ~two
    Cmono = Mc*c;
else
    nc = size(Y,2);
    Cmono = zeros((d+1)^2,nc);
    for col = 1:nc
        Cpq = reshape(c(:,col),d+1,d+1).';           % (p,q), rows were (p-1)(d+1)+q
        D = Mc*Cpq*Mc.';                             % (r,u)
        Cmono(:,col) = reshape(D.',[],1);            % (r-1)(d+1)+u
    end
end
end

function T = chebV(x,d)
% T(:,p+1) = T_p(x), p = 0..d
x = x(:);   T = zeros(numel(x),d+1);    T(:,1) = 1;
if d>=1,    T(:,2) = x;     end
for p = 2:d,    T(:,p+1) = 2*x.*T(:,p) - T(:,p-1);   end
end

function Phi = design2(T1,T2,d)
% rows (i,j) with i fastest, columns (p-1)(d+1)+q: Phi = T_p(x_i) T_q(y_j)
Phi = zeros(size(T1,1)*size(T2,1),(d+1)^2);
for q = 1:d+1
    for p = 1:d+1
        Phi(:,(p-1)*(d+1)+q) = reshape(T1(:,p)*T2(:,q)',[],1);
    end
end
end

function M = cheb2mono(d,alpha,beta)
% sum_p c_p T_p(alpha s + beta) = sum_r (M c)_r s^r; column p+1 of M holds T_p(alpha s + beta)
M = zeros(d+1,d+1);
Tprev = 1;  M(1,1) = 1;
if d>=1
    Tcur = [beta, alpha];   M(1:2,2) = Tcur(:);
    for p = 2:d
        Tnext = 2*conv(Tcur,[beta, alpha]);
        Tnext(1:numel(Tprev)) = Tnext(1:numel(Tprev)) - Tprev;
        M(1:numel(Tnext),p+1) = Tnext(:);
        Tprev = Tcur;   Tcur = Tnext;
    end
end
end

function [X,res,iter] = cg_sopvar_pd(A,B,tol,maxiter,X0,lambda,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,res,iter] = cg_sopvar_pd(A,B,tol,maxiter,X0,lambda,verbose) solves
% (A + lambda*I) X = B for a self-adjoint positive definite 'sopvar'
% operator A and a 'sopvar' B with a single input (scalar R input, i.e. B
% represents a function of the spatial variables), using the conjugate
% gradient method directly on the system. Unlike cg_sopvar, which applies CG
% to the normal equations A'*A X = A'*B for a general A, this needs no
% adjoint and only one application of A per iteration.
%
% INPUT
% - A:       sopvar, L2^p[S] -> L2^p[S], self-adjoint (A = A') and positive
%            definite. This is NOT checked: for a general A the iteration
%            can break down or diverge, use cg_sopvar instead;
% - B:       sopvar, R -> L2^p[S] (no input variables, dims(2)=1);
% - tol:     (optional) stop when ||B-A*X|| <= tol*||B||. Default 1e-8;
% - maxiter: (optional) maximum iterations. Default 10;
% - X0:      (optional) initial guess, sopvar R -> L2^p. Default 0.
% - lambda:  (optional) shift >= 0, solving (A+lambda*I) X = B, which
%            keeps the system positive definite and improves its
%            conditioning. Default 0;
% - verbose: (optional) print the residual norm of every iteration to the
%            command window. Default true.
%
% OUTPUT
% - X:       sopvar R -> L2^p[S], the solution;
% - res:     iterations+1 x 1 array of residual norms ||B - (A+lambda*I)*X_k||;
% - iter:    number of iterations performed.
%
% NOTES
% The inner product of two functions U,V (sopvar with a single input) is
% <U,V> = U'*V, an R -> R sopvar whose only parameter is the scalar
% integral of U(s)'V(s). All products use sopvar mtimes/ctranspose/plus/
% minus. Iterates grow in monomial degree with each application of A; no
% truncation is performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - cg_sopvar_pd
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
% Initial coding AT - 10/15/2026
% AT, 10/15/2026: Warning text also names rounding as a cause: near convergence
%                 <p,A*p> is computed with absolute error ~1e-10 (monomial
%                 inner products) and can come out <= 0 for a PD operator.

if nargin<3 || isempty(tol),     tol = 1e-8;     end
if nargin<4 || isempty(maxiter), maxiter = 10;   end
if nargin<6 || isempty(lambda),  lambda = 0;     end
if nargin<7 || isempty(verbose), verbose = true; end
if ~isscalar(lambda) || lambda<0,  error('lambda must be a scalar >= 0.');  end
if ~isa(A,'sopvar') || ~isa(B,'sopvar')
    error('A and B should be sopvar objects.')
end

% Standard CG for (A+lambda*I) X = B:
%   r = B - (A+lambda*I)*X       residual
%   p                            search direction
%   gam = <r,r>
if nargin<5 || isempty(X0)
    X = 0*B;   r = B;                   % zero initial guess
else
    X = X0;   r = B - A*X;
    if lambda>0,  r = r - lambda*X;  end
end
p = r;
gam = ip(r,r);
nB = sqrt(ip(B,B));
res = sqrt(gam);
iter = 0;
if verbose
    fprintf('cg_sopvar_pd: A is %d x %d, tol = %.1e, maxiter = %d\n',A.dims(1),A.dims(2),tol,maxiter);
    fprintf('%5s %14s %10s\n','iter','||B-A*X||','time (s)');
    fprintf('%5d %14.4e %10.2f\n',0,res(1),0);
    t0 = tic;
end

while iter<maxiter && sqrt(gam)>tol*nB && gam>0
    q = A*p;
    if lambda>0,  q = q + lambda*p;  end
    pq = ip(p,q);
    if pq<=0
%        warning('cg_sopvar_pd:notPD','<p,A*p> = %g <= 0 at iteration %d: A is not positive definite.',pq,iter+1);% AT, 10/15/2026 (was)
        warning('cg_sopvar_pd:notPD','<p,A*p> = %g <= 0 at iteration %d: A is not positive definite, or the residual reached rounding level (||r||/||B|| = %.1e).',pq,iter+1,sqrt(gam)/nB);% AT, 10/15/2026
        break
    end
    alpha = gam/pq;                     % alpha_k = <r_k,r_k>/<p_k,A p_k>
    X = X + alpha*p;
    r = r - alpha*q;                    % residual B - (A+lambda*I)*X
    gam_new = ip(r,r);
    p = r + (gam_new/gam)*p;            % conjugate search direction
    gam = gam_new;
    iter = iter+1;
    res(iter+1,1) = sqrt(gam);
    if verbose
        fprintf('%5d %14.4e %10.2f\n',iter,res(iter+1),toc(t0));
    end
end
if verbose
    if sqrt(gam)<=tol*nB
        fprintf('cg_sopvar_pd: converged in %d iterations.\n',iter);
    else
        fprintf('cg_sopvar_pd: stopped after %d iterations, ||r||/||B|| = %.2e.\n',iter,sqrt(gam)/nB);
    end
end
end


function val = ip(U,V)
% <U,V> = U'*V, an R -> R sopvar with a scalar parameter.
W = U'*V;
val = full(W.params{1});
end

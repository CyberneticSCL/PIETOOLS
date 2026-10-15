function [X,res,iter] = cg_copvar_pd(A,B,tol,maxiter,X0,lambda,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,res,iter] = cg_copvar_pd(A,B,tol,maxiter,X0,lambda,verbose) solves
% (A + lambda*I) X = B for a self-adjoint positive definite 'copvar'
% operator A and a 'copvar' B with a single input (an M x 1 container of
% sopvar blocks R -> L2, i.e. B represents a stacked vector of functions),
% using the conjugate gradient method directly on the system. Unlike
% cg_copvar, which applies CG to the normal equations A'*A X = A'*B for a
% general A, this needs no adjoint and only one application of A per
% iteration.
%
% INPUT
% - A:       copvar, M x M grid of sopvar blocks, a self-adjoint (A = A')
%            positive definite map of the stacked space L2^{p_1}[s^1] x ...
%            x L2^{p_M}[s^M] into itself. Neither property is checked: for a
%            general A the iteration can break down or diverge, use
%            cg_copvar instead;
% - B:       copvar, M x 1, blocks R -> L2^{p_i}[s^i] (no input variables);
% - tol:     (optional) stop when ||B-A*X|| <= tol*||B||. Default 1e-8;
% - maxiter: (optional) maximum iterations. Default 10;
% - X0:      (optional) initial guess, M x 1 copvar. Default 0;
% - lambda:  (optional) shift >= 0, solving (A+lambda*I) X = B, which
%            keeps the system positive definite and improves its
%            conditioning. Default 0;
% - verbose: (optional) print the residual norm of every iteration to the
%            command window. Default true.
%
% OUTPUT
% - X:       M x 1 copvar, the solution;
% - res:     iterations+1 x 1 array of residual norms ||B - (A+lambda*I)*X_k||;
% - iter:    number of iterations performed.
%
% NOTES
% The inner product of two stacked functions U,V (M x 1 copvars with a
% single input) is <U,V> = U'*V, a 1 x 1 copvar whose block is an R -> R
% sopvar with the scalar sum_i int U_i(s)'V_i(s); a structurally zero
% block means <U,V> = 0. All products use copvar mtimes/ctranspose/plus/
% minus. Iterates grow in monomial degree with each application of A; no
% truncation is performed.
%
% Same algorithm as cg_sopvar_pd, see there.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - cg_copvar_pd
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
if ~isa(A,'copvar') || ~isa(B,'copvar')
    error('A and B should be copvar objects.')
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
    fprintf('cg_copvar_pd: A is %d x %d blocks, tol = %.1e, maxiter = %d\n',size(A,1),size(A,2),tol,maxiter);
    fprintf('%5s %14s %10s\n','iter','||B-A*X||','time (s)');
    fprintf('%5d %14.4e %10.2f\n',0,res(1),0);
    t0 = tic;
end

while iter<maxiter && sqrt(gam)>tol*nB && gam>0
    q = A*p;
    if lambda>0,  q = q + lambda*p;  end
    pq = ip(p,q);
    if pq<=0
%        warning('cg_copvar_pd:notPD','<p,A*p> = %g <= 0 at iteration %d: A is not positive definite.',pq,iter+1);% AT, 10/15/2026 (was)
        warning('cg_copvar_pd:notPD','<p,A*p> = %g <= 0 at iteration %d: A is not positive definite, or the residual reached rounding level (||r||/||B|| = %.1e).',pq,iter+1,sqrt(gam)/nB);% AT, 10/15/2026
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
        fprintf('cg_copvar_pd: converged in %d iterations.\n',iter);
    else
        fprintf('cg_copvar_pd: stopped after %d iterations, ||r||/||B|| = %.2e.\n',iter,sqrt(gam)/nB);
    end
end
end


function val = ip(U,V)
% <U,V> = U'*V, a 1 x 1 copvar with an R -> R sopvar block (or a
% structurally zero block, for which the inner product is 0).
W = U'*V;
blk = W.C{1,1};
if isempty(blk),  val = 0;  else,  val = full(blk.params{1});  end
end

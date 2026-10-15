function [X,res,iter] = cg_copvar(A,B,tol,maxiter,X0,verbose,lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,res,iter] = cg_copvar(A,B,tol,maxiter,X0,verbose,lambda) solves
% A*X = B for a general 'copvar' operator A and a 'copvar' B with a single
% input (an M x 1 container of sopvar blocks R -> L2, i.e. B represents a
% stacked vector of functions), using conjugate gradient on the normal
% equations (A'*A + lambda*I) X = A'*B (CGLS with Tikhonov term).
%
% INPUT
% - A:       copvar, M x N grid of sopvar blocks, the map from the stacked
%            input spaces L2^{p_j}[s^j] to the stacked outputs L2^{q_i}[s^i];
% - B:       copvar, M x 1, blocks R -> L2^{q_i}[s^i] (no input variables);
% - tol:     (optional) stop when ||A'r|| <= tol*||A'B||. Default 1e-8;
% - maxiter: (optional) maximum iterations. Default 10;
% - X0:      (optional) initial guess, N x 1 copvar. Default 0;
% - verbose: (optional) print the residual and gradient norm of every
%            iteration to the command window. Default true;
% - lambda:  (optional) Tikhonov parameter >= 0, solving
%            min ||A*X-B||^2 + lambda*||X||^2. Default 0.
%
% OUTPUT
% - X:       N x 1 copvar, the (regularized least-squares) solution;
% - res:     iterations+1 x 1 array of residual norms ||B - A*X_k||;
% - iter:    number of iterations performed.
%
% NOTES
% The inner product of two stacked functions U,V (M x 1 copvars with a
% single input) is <U,V> = U'*V, a 1 x 1 copvar whose block is an R -> R
% sopvar with the scalar sum_i int U_i(s)'V_i(s); a structurally zero
% block means <U,V> = 0. The adjoint of A is ctranspose(A). All products
% use copvar mtimes/ctranspose/plus/minus. Iterates grow in monomial degree
% with each application of A; no truncation is performed.
%
% Same algorithm as cg_sopvar, see there for the regularization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - cg_copvar
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
% Initial coding AT - 10/14/2026
% Version of cg_sopvar for copvar containers. Unlike cg_sopvar, the initial
% X is 0*(A'*B) rather than 0*B, since X lives in the input space of A, which
% differs from the output space of B unless A is square with matching blocks.
% AT, 10/15/2026: Warnings when the iteration does not converge (maxiter reached
%                 with ||grad|| > tol*||A'B||), and when it breaks down: the
%                 monomial inner products <s,s>, ||r||^2 or <p,(A'A+lambda*I)p>
%                 come out <= 0 for large bases (Hilbert-type Gram matrix), in
%                 which case the result is not trustworthy.

if nargin<3 || isempty(tol),     tol = 1e-8;     end
if nargin<4 || isempty(maxiter), maxiter = 10;   end
if nargin<6 || isempty(verbose), verbose = true; end
if nargin<7 || isempty(lambda),  lambda = 0;     end
if ~isscalar(lambda) || lambda<0,  error('lambda must be a scalar >= 0.');  end
if ~isa(A,'copvar') || ~isa(B,'copvar')
    error('A and B should be copvar objects.')
end

At = A';                                % adjoint wrt the L2 inner product
if nargin<5 || isempty(X0)
    r = B;                              % zero initial guess
    s = At*r;
    X = 0*s;                            % X lives in the input space of A
else
    X = X0;   r = B - A*X;
    s = At*r;
end
% Gradient of the regularized functional: A'r - lambda*X.
if lambda>0,  s = s - lambda*X;  end
p = s;
gam = ip(s,s);
gam0 = ip(At*B,At*B);
res = sqrt(ip(r,r));
iter = 0;
if verbose
    fprintf('cg_copvar: A is %d x %d blocks, tol = %.1e, maxiter = %d\n',size(A,1),size(A,2),tol,maxiter);
    fprintf('%5s %14s %14s %10s\n','iter','||B-A*X||','||grad||','time (s)');
    fprintf('%5d %14.4e %14.4e %10.2f\n',0,res(1),sqrt(gam),0);
    t0 = tic;
end

% Description of a numerical breakdown, if any.                     % AT, 10/15/2026
bad = '';                                                                   % AT, 10/15/2026
while iter<maxiter && sqrt(gam)>tol*sqrt(gam0) && gam>0
    q = A*p;
    % <p,(A'A+lambda*I)p> = ||A*p||^2 + lambda*||p||^2.
    den = ip(q,q) + lambda*ip(p,p);                                         % AT, 10/15/2026
    % <p,(A'A+lambda*I)p> must be > 0.                                   % AT, 10/15/2026
    if den<=0                                                               % AT, 10/15/2026
        bad = sprintf('<p,(A''A+lambda*I)p> = %g <= 0 at iteration %d',den,iter+1);% AT, 10/15/2026
        break                                                               % AT, 10/15/2026
    end                                                                     % AT, 10/15/2026
    alpha = gam/den;                                                        % AT, 10/15/2026
    X = X + alpha*p;
    r = r - alpha*q;                    % residual B - A*X
    s = At*r;
    if lambda>0,  s = s - lambda*X;  end
    gam_new = ip(s,s);
    % <s,s> must be >= 0.                                                % AT, 10/15/2026
    if gam_new<0                                                            % AT, 10/15/2026
        bad = sprintf('<s,s> = %g < 0 at iteration %d',gam_new,iter+1);     % AT, 10/15/2026
        break                                                               % AT, 10/15/2026
    end                                                                     % AT, 10/15/2026
    p = s + (gam_new/gam)*p;            % conjugate search direction
    gam = gam_new;
    iter = iter+1;
    rr = ip(r,r);                                                           % AT, 10/15/2026
    % ||r||^2 must be >= 0.                                              % AT, 10/15/2026
    if rr<0 && isempty(bad)                                                 % AT, 10/15/2026
        bad = sprintf('||r||^2 = %g < 0 at iteration %d',rr,iter);          % AT, 10/15/2026
    end                                                                     % AT, 10/15/2026
    res(iter+1,1) = sqrt(max(rr,0));                                        % AT, 10/15/2026
    if verbose
        fprintf('%5d %14.4e %14.4e %10.2f\n',iter,res(iter+1),sqrt(gam),toc(t0));
    end
end
if ~isempty(bad)                                                            % AT, 10/15/2026
    warning('cg_copvar:breakdown','%s: inner products lost positivity (rounding in the monomial basis), result not trustworthy.',bad);% AT, 10/15/2026
elseif ~(sqrt(gam)<=tol*sqrt(gam0))                                         % AT, 10/15/2026
    warning('cg_copvar:noConvergence','Not converged after %d iterations: ||grad||/||A''B|| = %.2e > tol = %.1e.',iter,sqrt(gam/gam0),tol);% AT, 10/15/2026
end                                                                         % AT, 10/15/2026
if verbose
    if isempty(bad) && sqrt(gam)<=tol*sqrt(gam0)
        fprintf('cg_copvar: converged in %d iterations.\n',iter);
    else
        fprintf('cg_copvar: stopped at maxiter = %d, ||grad||/||A''B|| = %.2e.\n',maxiter,sqrt(gam/gam0));
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

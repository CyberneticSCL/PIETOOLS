function [X,res,iter] = cg_sopvar(A,B,tol,maxiter,X0,lambda,verbose)        % AT, 10/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,res,iter] = cg_sopvar(A,B,tol,maxiter,X0,verbose,lambda) solves A*X = B for a general
% 'sopvar' operator A and a 'sopvar' B with a single input 
% (scalar R input, i.e. B represents a function of the spatial variables), 
% using conjugate gradient on the normal equations 
% (A'*A + lambda I)*X = A'*B (CGLS).
%
% INPUT
% - A:       sopvar, L2^p[S1,S3] -> L2^q[S2,S3];
% - B:       sopvar, R -> L2^q[S2,S3] (no input variables, dims(2)=1);
% - tol:     (optional) stop when ||A'r|| <= tol*||A'B||. Default 1e-8;
% - maxiter: (optional) maximum iterations. Default 10;
% - X0:      (optional) initial guess, sopvar R -> L2^p. Default 0*B.
% - lambda:  (optional) Tikhonov parameter >= 0. Solves the regularized
%            problem min ||A*X-B||^2 + lambda*||X||^2, i.e.
%            (A'*A + lambda*I) X = A'*B. Default 0 (no regularization).
% - verbose: (optional) print the residual and gradient norm of every
%            iteration to the command window. Default true.
%
% OUTPUT
% - X:       sopvar R -> L2^p[S1,S3], (least-squares) solution;
% - res:     iterations+1 x 1 array of residual norms ||B - A*X_k||;
% - iter:    number of iterations performed.
%
% NOTES
% The inner product of two functions U,V (sopvar with a single input) is
% <U,V> = U'*V, an R -> R sopvar whose only parameter is the scalar
% integral of U(s)'V(s); the adjoint of A is then ctranspose(A). All
% products use sopvar mtimes/ctranspose/plus/minus. Iterates grow in
% monomial degree with each application of A; no truncation is performed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - cg_sopvar
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
% Initial coding AT - 10/13/2026

% AT, 10/15/2026: Warnings when the iteration does not converge (maxiter reached
%                 with ||grad|| > tol*||A'B||), and when it breaks down: the
%                 monomial inner products <s,s>, ||r||^2 or <p,(A'A+lambda*I)p>
%                 come out <= 0 for large bases (Hilbert-type Gram matrix), in
%                 which case the result is not trustworthy.

if nargin<3 || isempty(tol),     tol = 1e-8;     end
if nargin<4 || isempty(maxiter), maxiter = 10;  end               
if nargin<6 || isempty(lambda),   lambda = 0;      end       
if nargin<7 || isempty(verbose),  verbose = true;   end                     
if ~isscalar(lambda) || lambda<0,  error('lambda must be a scalar >= 0.');  end
if ~isa(A,'sopvar') || ~isa(B,'sopvar')
    error('A and B should be sopvar objects.')
end
% Applying standard CG method for \|A X - B\|^2 + lambda ||X||^2 or
% equvalently for (A'* A + lambda I)* X = A' * B
% A' * A is self-adjoint
% A' * B is a new residual
%
% % % % % % % % % % % % % variables in the code % % % % % % % % % % % % % 
% X is a solution 
% r is a residual A*X - B
% s is a A'*A*X - A'*B
% p is a gradient direction 
% gam is a scalar of s'*s
% % % % % % % % % % % % % % % ALGORITHM % % % % % % % % % % % % % % % % %

if nargin<5 || isempty(X0)
    X = 0*(A'*B);   r = B;                    % zero initial guess
else
    X = X0;   r = B - A*X;
end

At = A';                                % adjoint wrt the L2 inner product
s = At*r;

% Gradient of the regularized functional: A'r - lambda*X.               
if lambda>0,  s = s - lambda*X;  end                                      
p = s;
gam = ip(s,s);
gam0 = ip(At*B,At*B);
res = sqrt(ip(r,r));
iter = 0;
if verbose
    fprintf('cg_sopvar: A is %d x %d, tol = %.1e, maxiter = %d\n',A.dims(1),A.dims(2),tol,maxiter);
    fprintf('%5s %14s %14s %10s\n','iter','||B-A*X||','||A''(B-A*X)||','time (s)');
    fprintf('%5d %14.4e %14.4e %10.2f\n',0,res(1),sqrt(gam),0);             
    t0 = tic;                                                               
end                                                                         

fprintf("")
% Description of a numerical breakdown, if any.                     % AT, 10/15/2026
bad = '';                                                                   % AT, 10/15/2026
while iter<maxiter && sqrt(gam)>tol*sqrt(gam0) && gam>0
    q = A*p;
    %   alpha = gam/ip(q,q); % alpha_k = r_k^T r_k/ (p_k, A'*A p_k)         
    % <p,(A'A+lambda*I)p> = ||A*p||^2 + lambda*||p||^2.                 
    den = ip(q,q) + lambda*ip(p,p);                                         % AT, 10/15/2026
    % <p,(A'A+lambda*I)p> must be > 0.                                   % AT, 10/15/2026
    if den<=0                                                               % AT, 10/15/2026
        bad = sprintf('<p,(A''A+lambda*I)p> = %g <= 0 at iteration %d',den,iter+1);% AT, 10/15/2026
        break                                                               % AT, 10/15/2026
    end                                                                     % AT, 10/15/2026
    alpha = gam/den;                                                        % AT, 10/15/2026
    X = X + alpha*p;  % X_k+1 = X_k + alpha p_k
    r = r - alpha*q;  % residual B - A X_{k+1}
    s = At*r; % residual A' * B - A' * A * X_{k+1}
    % Regularized gradient A'r - lambda*X.                               
    if lambda>0,  s = s - lambda*X;  end                                   
    gam_new = ip(s,s); % beta_k = ||r_{k+1}||^2/||r_k||^2
    % <s,s> must be >= 0.                                                % AT, 10/15/2026
    if gam_new<0                                                            % AT, 10/15/2026
        bad = sprintf('<s,s> = %g < 0 at iteration %d',gam_new,iter+1);     % AT, 10/15/2026
        break                                                               % AT, 10/15/2026
    end                                                                     % AT, 10/15/2026
    p = s + (gam_new/gam)*p;  % conjugate search direction
    % p_{k+1} = s_{k+1} - beta_k p_k
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
    warning('cg_sopvar:breakdown','%s: inner products lost positivity (rounding in the monomial basis), result not trustworthy.',bad);% AT, 10/15/2026
elseif ~(sqrt(gam)<=tol*sqrt(gam0))                                         % AT, 10/15/2026
    warning('cg_sopvar:noConvergence','Not converged after %d iterations: ||A''r||/||A''B|| = %.2e > tol = %.1e.',iter,sqrt(gam/gam0),tol);% AT, 10/15/2026
end                                                                         % AT, 10/15/2026
if verbose                                                                
    if isempty(bad) && sqrt(gam)<=tol*sqrt(gam0)                           
        fprintf('cg_sopvar: converged in %d iterations.\n',iter);          
    else                                                                   
        fprintf('cg_sopvar: stopped at maxiter = %d, ||A''r||/||A''B|| = %.2e.\n',maxiter,sqrt(gam/gam0));% AT, 10/13/2026
    end                                                                     
end                                                                         
end


function val = ip(U,V)
% <U,V> = U'*V, an R -> R sopvar with a scalar parameter.
W = U'*V;
val = full(W.params{1});
end

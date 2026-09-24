function [F,J,q] = bm_resid(w,P,rv,skipF)                                   % CC, 09/22/2026
% skipF (optional): true returns F = [] and q = [] and computes ONLY J.
% WHY.  bm_lm2 accepts a step, sets F = Fn (already evaluated in the damping
% loop), and then calls [F,J] = bm_resid(w,...) purely to get J at the new
% point -- recomputing an F it already holds.  MEASURED time split of one LM
% iteration (1-D heavy, rank 2/block): (F,J) recompute 41.6%, trial F 32.9%,
% chol+solve 20.1%, forming J'J 5.4% -- so F alone is ~17% of the step and the
% residual/Jacobian evaluations are 74.5% of it, against 20% for the linear
% algebra.  J needs only Y, never q, so skipping F also skips the whole
% Y*Y' -> q assembly and the Ssym/W products, not just the final norm.
% Preconditioned residual  F = W*(A(YY') - b)  and its analytic Jacobian.
% With the preconditioner of bm_setup, ||F|| is the Frobenius distance from
% Q = blkdiag(Y_i Y_i') to the affine set {A(X) = b}.
% w = [ z ; Y_1(:) ; ... ; Y_B(:) ],  Y_i is Ns(i) x rv(i)  (rv(i)=0 -> Q_i=0).
if nargin<4||isempty(skipF), skipF = false; end                             % CC, 09/22/2026
B = numel(P.Ns);
q = zeros(P.Ntot,1);
k = P.Kf;
% w still stores the free coordinates FIRST; only their destination in q      % CC, 09/23/2026
% changes, because a program built with lpivar interleaves them between the   % CC, 09/23/2026
% Gram blocks instead of placing them in a prefix.  P.free is 1:Kf whenever    % CC, 09/23/2026
% the old contiguous layout holds, so this is a no-op there.                  % CC, 09/23/2026
%if P.Kf>0, q(1:P.Kf) = w(1:P.Kf); end                                       % CC, 09/23/2026 (was)
if P.Kf>0, q(P.free) = w(1:P.Kf); end                                       % CC, 09/23/2026
Y = cell(1,B);
for i=1:B
    N = P.Ns(i); r = rv(i);
    if r==0, Y{i} = zeros(N,0); continue; end
    Y{i} = reshape(w(k+(1:N*r)),N,r);  k = k+N*r;
    if skipF, continue; end             % J needs Y only, never q            % CC, 09/22/2026
    Qi = Y{i}*Y{i}';
    q(P.rows{i}) = Qi(:);
end
if skipF                                                                    % CC, 09/22/2026
    F = [];  q = [];                                                        % CC, 09/22/2026
else                                                                        % CC, 09/22/2026
    F = P.W*(P.Ssym*q - P.bf);
end                                                                         % CC, 09/22/2026
if nargout<2 && ~skipF, return; end                                         % CC, 09/22/2026
Jc = cell(1,B+1);
%if P.Kf>0, Jc{1} = P.W*P.Ssym(:,1:P.Kf); else ...                          % CC, 09/23/2026 (was)
if P.Kf>0, Jc{1} = P.W*P.Ssym(:,P.free); else, Jc{1} = zeros(P.mres,0); end % CC, 09/23/2026
for i=1:B
    N = P.Ns(i); r = rv(i);
    if r==0, Jc{i+1} = zeros(P.mres,0); continue; end
    Jc{i+1} = 2*(P.W*(P.Ssym(:,P.rows{i})*kron(sparse(Y{i}),speye(N))));
end
J = [Jc{:}];
end

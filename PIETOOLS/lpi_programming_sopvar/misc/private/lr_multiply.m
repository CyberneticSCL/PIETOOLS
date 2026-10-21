function [Aout,Bout] = lr_multiply(L,A,B,R)                                 % MMP, 09/30/2026
% function [Aout,Bout] = lr_multiply(L,A,B,R,K)                             % MMP, 09/30/2026 (was)
% function [Aout,Bout] = lr_multiply(L,A,B,R)                               % MMP, 09/26/2026 (was)
% Given vec(C(d)) = A + B'*d, return the coefficients of vec(L*C(d)*R), using
% vec(L*C*R) = (R' kron L)*vec(C).
%
% The constant term is handled by reshaping, which avoids the Kronecker
% (was) product entirely. For B the product is formed explicitly, as in     % MMP, 09/30/2026 (was)
% (was) @sdopvar/plus: here L and R are monomial selection matrices with a single % MMP, 09/30/2026 (was)
% product entirely. For B the product is formed explicitly:                 % MMP, 09/30/2026
% here L and R are monomial selection matrices with a single                % MMP, 09/30/2026
% nonzero per row, so the Kronecker product has only nnz(L)*nnz(R) nonzeros
% and applying it to all decision variables at once is far cheaper than
% looping over the rows of B, of which there may be many thousands.
%
% (was) Optional K = kron(R.',L).', for a caller applying the same L and R to     % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% (was) several coefficient sets (one per gamma cell): forming it once per pair   % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% (was) instead of once per cell removed ~3 of 11.7 s of a 2-D poscopvar build.   % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
%
% MMP, 09/26/2026: optional fifth input K, as above. Without it the
%                  behaviour is unchanged.
% MMP, 09/26/2026: Skip B when only Aout is requested: 'unpack_sheets' now
%                  maps B itself and calls this for A alone, so forming K
%                  there would undo the saving. No caller passes K any more
%                  ('copquadvar' did, above); the input is kept.
% MMP, 09/30/2026: Fifth input K removed. No caller passed it (previous
%                  entry), so the first 09/26/2026 entry no longer holds;
%                  Bout is B*kron(R.',L).' again, the product the K default
%                  formed. The removed default ran only on the unfused path
%                  of 'unpack_sheets' and in 'sopquadvar', both 0 calls in
%                  the w1-w3 builds. Help: no longer "as in @sdopvar/plus",
%                  which maps its bases with 'apply_basis_map' (via
%                  plus_batch) since 09/07/2026.

X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;
Aout = Y(:);
if nargout<2,  return;  end     % B unused: no Kronecker product            % MMP, 09/26/2026

% Bout = B*kron(R.',L).';                                                   % MMP, 09/26/2026 (was)
% if nargin<5 || isempty(K),  K = kron(R.',L).';  end                         % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
% Bout = B*K;                                                                 % MMP, 09/26/2026 % MMP, 09/30/2026 (was)
Bout = B*kron(R.',L).';     % the product the K default formed              % MMP, 09/30/2026

end

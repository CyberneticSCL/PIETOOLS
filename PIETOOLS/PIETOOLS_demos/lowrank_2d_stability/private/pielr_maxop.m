function [m,npar] = pielr_maxop(X)                                          % CC, 09/23/2026
% PIELR_MAXOP  max|coefficient| over EVERY parameter of a substituted
% container operator ('opvar'), in any dimension, and NPAR the number of
% parameters visited.
%
% Every block and every parameter cell is visited, because a silently missed
% cell would let a nonzero residual read as zero -- the one failure mode that
% turns a wrong answer into a passing one.
%
% MMP, 11/26/2026: the branch for the opvar, dopvar, opvar2d and dopvar2d
%   classes deleted, with its help (see the end of the code).

% BEGIN CC, 11/08/2026: a 'copvar' (substituted container): the largest
% |coefficient| over every parameter of every 'sopvar' block of X.C.  A block
% holding [] is structurally zero and has no parameters.  No enumeration
% assert: the block cell and the parameter cells are walked whole, so a
% block or multi-index cannot be skipped the way a hard-coded field list can.
if isa(X,'opvar')
    m = 0;  npar = 0;
    for i = 1:numel(X.C)
        B = X.C{i};
        if isempty(B),  continue,  end
        for k = 1:numel(B.params)
            v = B.params{k};
            if ~isempty(v),  m = max(m,full(max(abs(v(:)))));  end
            npar = npar + 1;
        end
    end
    return
end
% END CC, 11/08/2026
% BEGIN MMP, 11/26/2026: deleted, CC 09/23/2026: the enumeration of the
% properties of opvar and dopvar (4 properties, 6 leaves) and of opvar2d and
% dopvar2d (16, 36), its leaf-count assert, and the local function mxc.
% Those classes left the path (f83486c8); once dopvar is renamed dopvar, a
% dopvar would reach that enumeration and fail the assert. Any other class
% is an input error.
error('pielr_maxop:input','expected an opvar (a substituted container), got %s',class(X));
% END MMP, 11/26/2026
end

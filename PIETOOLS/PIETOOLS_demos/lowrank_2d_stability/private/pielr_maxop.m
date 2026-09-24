function [m,npar] = pielr_maxop(X)                                          % CC, 09/23/2026
% PIELR_MAXOP  max|coefficient| over EVERY parameter cell of a PI operator,
% in either dimension, with the leaf count ASSERTED.
%
% Generalises opcheck_2d's local mxop2d to 1-D.  opvar has 4 parameter
% properties and 6 leaves (P, Q1, Q2, R.R0, R.R1, R.R2); opvar2d has 16
% properties and 36 leaves (R00,R0x,R0y,R02, Rx0,Rxx{3},Rxy,Rx2{3}, Ry0,Ryx,
% Ryy{3},Ry2{3}, R20,R2x{3},R2y{3},R22{3,3}).
%
% The enumeration is taken from properties() rather than hard-coded, and the
% leaf count is asserted, because a silently missed cell would let a nonzero
% residual read as zero -- the one failure mode that turns a wrong answer into
% a passing one.  (opcheck_2d's own reasoning; kept, and extended to 1-D,
% where the previous 1-D checker hard-coded the six names with no assert.)
%
% The 1-D R field is a STRUCT, not a cell, so it needs its own branch: an
% opvar's properties() returns 'R' once, and its three sub-fields are leaves.

if isa(X,'opvar') || isa(X,'dopvar')
    nprop_want = 4;   nleaf_want = 6;    skip = {'I','var1','var2','dim','dimdependent'};
elseif isa(X,'opvar2d') || isa(X,'dopvar2d')
    nprop_want = 16;  nleaf_want = 36;   skip = {'I','var1','var2','dim','dimdependent'};
else
    error('pielr_maxop:input','expected an opvar/opvar2d (or dopvar variant), got %s',class(X));
end
pr   = properties(X);
pars = pr(~ismember(pr,skip));
m = 0;  npar = 0;
for i = 1:numel(pars)
    v = X.(pars{i});
    if iscell(v)
        for k = 1:numel(v)
            m = max(m,mxc(v{k}));  npar = npar+1;
        end
    elseif isstruct(v)
        % the 1-D R = struct('R0',..,'R1',..,'R2',..)
        fn = fieldnames(v);
        for k = 1:numel(fn)
            m = max(m,mxc(v.(fn{k})));  npar = npar+1;
        end
    else
        m = max(m,mxc(v));  npar = npar+1;
    end
end
assert(numel(pars)==nprop_want && npar==nleaf_want, ...
    'pielr_maxop: %s enumeration changed (%d fields / %d leaves, expected %d / %d)', ...
    class(X),numel(pars),npar,nprop_want,nleaf_want);
end

% =========================================================================
function m = mxc(X)
if isempty(X), m = 0; return; end
if     isa(X,'double'),     C = X;
elseif isa(X,'polynomial'), C = X.coefficient;
elseif isa(X,'dpvar'),      C = X.C;
else,  m = NaN; return;                       % unknown leaf type: poison, do not hide
end
if isempty(C), m = 0; else, m = full(max(abs(C(:)))); end
if isempty(m), m = 0; end
end

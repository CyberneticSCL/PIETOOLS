function C = plus_batch(varargin)
% C = PLUS_BATCH(P1,P2,...,PN) is P1 + P2 + ... + PN for 'copvar'
% containers on one block grid. Fixed containers carry no decision variable
% lists, so the chain has nothing to reconcile and is used as is. A call
% with any 'cdopvar' operand dispatches to '@cdopvar/plus_batch'
% (InferiorClasses), which reconciles the lists once.
%
% Initial coding MMP, 10/01/2026

if ~all(cellfun(@(x) isa(x,'copvar'),varargin))
    error('plus_batch:badInput','Summands must be copvar or cdopvar objects.')
end
C = varargin{1};
for k = 2:nargin,   C = C + varargin{k};    end
end

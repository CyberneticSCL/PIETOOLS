function [r,res] = mosekopt(cmd,prob,varargin) %#ok<INUSL>
% [R,RES] = MOSEKOPT(CMD,PROB) - BENCHMARK STUB, not the MOSEK solver.
% Returns at once with a zero solution of the sizes PROB implies and
% status UNKNOWN, so that a run measures everything before the SDP solver
% (the parsing time). Put this folder first on the path of a parse-only
% benchmark session only.
% INPUT
% - cmd:  char, ignored;
% - prob: 'struct' of SEDUMI2MOSEK: fields c, a, bardim are read for sizes.
% OUTPUT
% - r:    0;
% - res:  'struct' with res.sol.itr.{xx,y,barx,bars} zero, prosta and
%         solsta 'UNKNOWN', res.info.MSK_DINF_OPTIMIZER_TIME = 0.

% Initial coding MMP, 11/06/2026: parse-only timing for the scaling plots.
n = 0;  if isfield(prob,'c'),       n = numel(prob.c);      end
m = 0;  if isfield(prob,'a'),       m = size(prob.a,1);     end
nb = 0; if isfield(prob,'bardim'),  nb = sum(prob.bardim.*(prob.bardim+1)/2);  end
res.sol.itr.xx = zeros(n,1);    res.sol.itr.y = zeros(m,1);
res.sol.itr.barx = zeros(nb,1); res.sol.itr.bars = zeros(nb,1);
res.sol.itr.prosta = 'UNKNOWN'; res.sol.itr.solsta = 'UNKNOWN';
res.info.MSK_DINF_OPTIMIZER_TIME = 0;
res.info.MSK_IINF_INTPNT_ITER = 0;
res.info.MSK_DINF_INTPNT_OPT_STATUS = 0;
r = 0;
end

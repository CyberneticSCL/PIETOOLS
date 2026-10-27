function [prog0,out] = cx_stock_capture(exec,PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG0,OUT] = CX_STOCK_CAPTURE(EXEC,PIE,ST) runs the UNMODIFIED stock
% executive 'PIETOOLS_<EXEC>' on PIE with settings ST, but with lpisolve
% shadowed by the capture stub of PIETOOLS_demos/cuadmm (cuadmm_shadow), and
% returns the program the executive handed to the solver, UNSOLVED.
%
% Used to pose the stock LPI at a FIXED gamma (cx_stock_pinned) without
% hand-transcribing the executive: the stock 1-D executives hard-declare
% gamma as an objective, and a transcription could not be checked against
% the original except by its numbers.
%
% OUT holds the executive's outputs, computed from the stub's zero solution:
% meaningless as values, returned only so the call does not error.
%
% The shadow is removed on exit, error or not (onCleanup). An executive that
% calls lpisolve more than once leaves its LAST program in PROG0.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global CENSUS_PROG
CENSUS_PROG = [];
w1 = warning('off','cuadmm_shadow:stubCalled');
w2 = warning('off','cuadmm_shadow:active');
cuadmm_shadow('lpisolve');
cleanup = onCleanup(@() restore(w1,w2));                                     %#ok<NASGU>
fn = ['PIETOOLS_' exec];
nout = abs(nargout(fn));
out = cell(1,nout);
evalc('[out{1:nout}] = feval(fn,PIE,st);');
prog0 = CENSUS_PROG;
if isempty(prog0)
    error('cx_stock_capture:none','%s never called lpisolve.',fn)
end
end

function restore(w1,w2)
cuadmm_shadow('off');
rehash;
warning(w1);    warning(w2);
end

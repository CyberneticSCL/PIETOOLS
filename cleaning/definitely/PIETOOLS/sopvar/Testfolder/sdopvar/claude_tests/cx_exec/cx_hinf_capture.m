function [prog0,perr] = cx_hinf_capture(exec,PIE,st,gain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG0,PERR] = CX_HINF_CAPTURE(EXEC,PIE,ST,GAIN) is cx_stock_capture with
% the executive's own fixed-gain argument: it runs the UNMODIFIED
% 'PIETOOLS_<EXEC>(PIE,ST,GAIN)' with lpisolve shadowed by the capture stub
% (cuadmm_shadow) and returns the program handed to the solver, unsolved.
%
% Why: the 2-D gain executives take a gain as third input and then build the
% LPI at that FIXED gamma (e.g. PIETOOLS_Hinf_gain_2D.m:182-186), which is
% exactly what the container transcription builds; cx_stock_capture calls
% with two inputs, i.e. the objective branch (gam a decision variable).
% Both are reported. 2-D programs are never solved here.
%
% The executive may error AFTER the capture while post-processing the
% stub's zero solution (e.g. lpigetsol on a numeric gam); PROG0 is still the
% program it built, and PERR carries the message ('' if none).
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
nout = abs(nargout(fn));    out = cell(1,nout);                             %#ok<NASGU>
perr = '';
try
    evalc('[out{1:nout}] = feval(fn,PIE,st,gain);');
catch ME
    perr = sprintf('[%s] %s',ME.identifier,ME.message);
end
prog0 = CENSUS_PROG;
if isempty(prog0)
    error('cx_hinf_capture:none','%s never called lpisolve (%s).',fn,perr)
end
end

function restore(w1,w2)
cuadmm_shadow('off');
rehash;
warning(w1);    warning(w2);
end

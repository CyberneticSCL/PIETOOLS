function prog_sol = lpisolve(prog,opts)
% LPISOLVE size gate of the library pass (baselines/2026-11-17_tiers; MMP,
% 11/17/2026): the repository LPISOLVE body behind two appdata caps,
% LB_NDV_CAP on the decision variables (rows of decvartable) and LB_M_CAP on
% the equality rows (columns of expr.At over the 'eq' expressions, the
% count of sossolve's 'Size:' line). Above a cap the program is not solved:
% the error 'LB_SIZE ndv=.. m=..' becomes the row's error note in lb_run,
% so the size is recorded. Without the appdata entries nothing is capped.
% lb_run puts this folder first on the path after lb_setup.
ndv = size(prog.decvartable,1);
iseq = strcmp(prog.expr.type,'eq');
m = sum(cellfun(@(A) size(A,2),prog.expr.At(iseq)));
cap = Inf;      if isappdata(0,'LB_NDV_CAP'),   cap = getappdata(0,'LB_NDV_CAP');    end
mcap = Inf;     if isappdata(0,'LB_M_CAP'),     mcap = getappdata(0,'LB_M_CAP');     end
if ndv > cap || m > mcap
    error('lb:size','LB_SIZE ndv=%d m=%d above the gate (ndv %g, m %g): not solved',ndv,m,cap,mcap);
end
% The repository body (lpi_programming/lpisolve.m, MMP 11/12/2026).
if nargin<2 || isempty(opts)
    opts = [];
    if isfield(prog,'settings') && isfield(prog.settings,'sos_opts'),   opts = prog.settings.sos_opts;  end
elseif isfield(opts,'sos_opts')
    opts = opts.sos_opts;
end
if isempty(opts),   prog_sol = sossolve(prog);
else,               prog_sol = sossolve(prog,opts);
end
end

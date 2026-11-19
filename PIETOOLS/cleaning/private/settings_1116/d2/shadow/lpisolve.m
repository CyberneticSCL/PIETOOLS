function prog_sol = lpisolve(prog,opts)
% LPISOLVE scratch shadow (sol1117/size, 11/17/2026): the repository
% LPISOLVE body behind a size gate. Prints '@@ LB_SIZE ndv=.. m=.. ineq=..':
% ndv = rows of decvartable (the first entry of sossolve's 'Size:' line on
% the LPI path, where every Gram entry is a declared decision variable);
% m = columns of expr.At over the 'eq' expressions (Atf = [expr.At{:}];
% processvars keeps the column count); ineq = expressions of any other type
% (0 on the LPI path; addextrasosvar would add variables otherwise).
% Appdata gates: LB_SIZE_ONLY true (never solve), LB_NDV_CAP (default 3e6),
% LB_M_CAP (default Inf); above a cap the error 'LB_SIZE ...' records the
% size instead of a solve. The 11/17 sizing of D2 reproduced the Sol
% 'Size:' line exactly (6,813,852 x 163,847).
ndv = size(prog.decvartable,1);
iseq = strcmp(prog.expr.type,'eq');
m = sum(cellfun(@(A) size(A,2),prog.expr.At(iseq)));
nin = sum(~iseq);
fprintf('@@ LB_SIZE ndv=%d m=%d ineq=%d\n',ndv,m,nin);
cap = 3e6;      if isappdata(0,'LB_NDV_CAP'),   cap = getappdata(0,'LB_NDV_CAP');    end
mcap = Inf;     if isappdata(0,'LB_M_CAP'),     mcap = getappdata(0,'LB_M_CAP');     end
only = isappdata(0,'LB_SIZE_ONLY') && getappdata(0,'LB_SIZE_ONLY');
if only || ndv > cap || m > mcap
    error('lb:size','LB_SIZE ndv=%d m=%d ineq=%d (caps ndv %g, m %g, size only %d)',ndv,m,nin,cap,mcap,only);
end
% The repository body (lpi_programming_sopvar/lpisolve.m, MMP 11/12/2026).
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

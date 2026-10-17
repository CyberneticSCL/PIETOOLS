function [prog,aux,hist] = lpi_solve_loop_sop(build,st,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,AUX,HIST] = LPI_SOLVE_LOOP_SOP(BUILD,ST[,LABEL]) the degree loop of
% the container executives. BUILD is a function handle
%
%   [prog,aux] = build(dw)
%
% that assembles the LPI with every slack weight raised by DW over the
% executive's rule and returns the program and a struct AUX (the executive's
% operators and the objective in AUX.obj, read by LPI_CLASSIFY_SOP). The
% loop solves with ST.sos_opts, classifies the exit, and while the status is
% 'unknown' (MOSEK status UNKNOWN: the signature of a slack degree that
% admits no certificate at any finite value, sopvar_lift_notes Sec. 10) or
% 'infeasible' (a slack too small for the operator reads as primal
% infeasibility when the objective is bounded below; for a feasibility LPI
% the raises cost solves and change nothing) and fewer than ST.sop.max_raise
% raises have been made, raises DW by one and
% rebuilds. HIST records every pass (dw, status, objective, numerr, pinf,
% dinf, decision count, rows, build and solve times); AUX is that of the
% last pass with AUX.classify and AUX.dw added.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3,    label = '';     end
st = exec_sop_settings(st);
hist = struct('dw',{},'status',{},'divergent',{},'obj',{},'numerr',{},'pinf',{},'dinf',{}, ...
              'ndv',{},'m',{},'t_build',{},'t_solve',{});
dw = 0;
for it = 0:st.sop.max_raise
    t0 = tic;   [prog,aux] = build(dw);     tb = toc(t0);
    t0 = tic;   prog = lpisolve(prog,st.sos_opts);  ts = toc(t0);
    obj = [];   if isfield(aux,'obj'),  obj = aux.obj;    end
    cl = lpi_classify_sop(prog,obj);
    m = 0;
    if isfield(prog,'expr') && isfield(prog.expr,'At')
        for e = 1:numel(prog.expr.At),  m = m + size(prog.expr.At{e},2);    end
    end
    hist(end+1) = struct('dw',dw,'status',cl.status,'divergent',cl.divergent,'obj',cl.obj, ...
                         'numerr',cl.numerr,'pinf',cl.pinf,'dinf',cl.dinf, ...
                         'ndv',numel(prog.decvartable),'m',m,'t_build',tb,'t_solve',ts); %#ok<AGROW>
    if st.sop.verbose
        fprintf('  [%s] dw %d: %s%s, objective %s, %d decision variables, %d rows, build %.1f s, solve %.1f s\n', ...
                label,dw,cl.status,tern(cl.divergent,' (divergent)',''),num2str(cl.obj,'%.8g'), ...
                numel(prog.decvartable),m,tb,ts);
    end
%   if ~strcmp(cl.status,'unknown') || it==st.sop.max_raise,   break,  end % MMP, 10/08/2026 (was)
    if ~any(strcmp(cl.status,{'unknown','infeasible'})) || it==st.sop.max_raise,  break,  end % MMP, 10/08/2026
    dw = dw+1;
end
aux.classify = cl;  aux.dw = dw;    aux.hist = hist;
end


function s = tern(c,a,b)
if c,   s = a;  else,   s = b;  end
end

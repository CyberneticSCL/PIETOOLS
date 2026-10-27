function R = cx_h2_struct2d(cases2d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = CX_H2_STRUCT2D(CASES2D) compares, WITHOUT SOLVING, the 2-D H2 programs
% built by the stock executive and by its container transcription:
%
%   stock      cx_stock_capture: the unmodified executive with lpisolve
%              shadowed, so the program is captured before any solve;
%   container  cx_<exec>(PIE,st,gam) at the case's fixed gam.
%
% Both shapes are cx_shape of the UNSOLVED program, so neither includes the
% slack blocks sossolve adds for 'ineq' rows; the comparison is like with
% like. Prints one row per case: shapes and assembly wall times. The stock
% program carries gam as a free objective variable (L130-134), the
% container a number, so Kf and ndv differ by one by construction.
%
% CASES2D: struct array with fields id, exec, setname, solver, plant, gam
% (the second output of cx_cases_h2).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
R = cell(1,numel(cases2d));
for ic = 1:numel(cases2d)
    c = cases2d(ic);
    st = cx_settings(c.setname,c.solver);
    PIE = cx_plant(c.plant{:});
    r = struct('id',c.id,'exec',c.exec,'err','');
    try
        t = tic;    prog0 = cx_stock_capture(c.exec,PIE,st);    r.tstock = toc(t);
        r.stock = cx_shape(prog0);
        t = tic;    prog = feval(['cx_' c.exec],PIE,st,c.gam);  r.tcx = toc(t);
        r.cx = cx_shape(prog);
        S = r.stock;    C = r.cx;
        fprintf(['  %-12s (structure only, gam %g) ndv %d/%d  Kf %d/%d  m %d/%d  Ks %s / %s' ...
                 '  | assembly stock %.1f s, cx %.1f s\n'],c.id,c.gam,S.ndv,C.ndv,S.Kf,C.Kf, ...
                 S.m,C.m,mat2str(S.Ks),mat2str(C.Ks),r.tstock,r.tcx);
    catch ME
        r.err = sprintf('[%s] %s (%s:%d)',ME.identifier,ME.message,ME.stack(1).name,ME.stack(1).line);
        fprintf('  %-12s ERROR %s\n',c.id,r.err);
    end
    R{ic} = r;
end
end

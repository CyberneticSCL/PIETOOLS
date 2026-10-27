function r = cx_stock_pinned(prog0,gam,sopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = CX_STOCK_PINNED(PROG0,GAM,SOPTS) poses a captured stock objective
% program (cx_stock_capture) as the FEASIBILITY question at gamma = GAM and
% solves it with cx_solve.
%
% The objective of every stock Hinf / H2 executive is one unit entry on the
% decision variable gamma. Pinning adds ONE equality, gamma - GAM = 0, and
% zeroes the objective; everything else is the program the executive built,
% solved through the same lpisolve/sossolve path (b normalized after the
% row is added, which the SDP-level bl_fix has to arrange by hand).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = find(prog0.objective);
if numel(j)~=1 || abs(full(prog0.objective(j))-1)>1e-12
    error('cx_stock_pinned:obj',['Expected one unit objective entry, found %d; '...
          'this program is not a min-gamma LPI.'],numel(j))
end
gname = prog0.decvartable{j};
gv = dpvar(gname);
prog = soseq(prog0,gv - gam);
prog.objective(:) = 0;
r = cx_solve(prog,sopts);
r.gam = gam;    r.gname = gname;
end

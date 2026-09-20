% DEMO_LOWRANK_2D  Heat-equation showcase: 2-D stability where SOLVING costs
% dramatically less than SETUP.
%
% This is a SCRIPT: everything it computes stays in your workspace afterwards,
% in particular the certificates certs.n1, certs.n2, ... (and cert1, rows),
% so you can run further tests on them -- e.g.
%     certN = pielr_certify(nb_rd2d(4), pielr_tensor(certs.n1,4));
%
% CC, 09/19/2026: converted from a function to a SCRIPT (maintainer request)
%                 so results persist in the workspace; the helpers it calls
%                 (nb_rd2d, set2d_deg, build_stab_2d_st2) moved from private/
%                 to the package folder, since scripts cannot see private/.
%
% System: n-state 2-D reaction-diffusion (heat) equation on the unit square,
%   x_t = x_{s1s1} + x_{s2s2} + 2*x, Dirichlet on all four edges (stable:
%   lam = 2 < lam* = 2*pi^2).  The n-state member is n identical copies, which
%   is exactly what makes the state ladder cheap: the certified n=1 face
%   replicates per monomial group (pielr_tensor) and certifies every larger n
%   with at most a tiny SDP -- measured, the full stock solve of the SAME n=1
%   program took 1707.3 s, and at n=4 the stock route needs a 66.8 GiB dense
%   normal matrix and is memory-infeasible on a 63.7 GiB machine.
%
% Every row is verified LIVE at the operator level (residual over all 36
% opvar2d cells < 1e-6 AND every Gram block PSD).  The script aborts loudly if
% any row fails.

% ---------------- config ----------------------------------------------------
nlist     = [1 2 3];   % state counts to run.  n = 4,5,6 also certify (measured)
                       % but the LPI ASSEMBLY needs roughly 8, 10, 13+ GiB of
                       % RAM and 6-10 min per row; extend when you have both.
DO_SEARCH = false;     % true : discover the n=1 face live (BM, minutes).
                       % false: use the shipped face_n1.mat and RE-VERIFY it
                       %        live against a freshly built n=1 program.
RUN_STOCK = false;     % true : also run the stock full solve at n=1 (~28 min
                       %        measured) for a live reference point.
% -----------------------------------------------------------------------------

here = fileparts(mfilename('fullpath'));
cd(here);              % 2-D conversion measured to fail from a cluttered cwd
pielr_path();

st = set2d_deg(4,[]);  % the demo's basis knob: eq degrees = LF degrees + 4

% ---- the n=1 face -----------------------------------------------------------
if DO_SEARCH
    fprintf('\n=== n = 1: cold discovery (DO_SEARCH = true) ===\n'); %#ok<UNRCH> % config constant
    PIE1 = nb_rd2d(1);
    o1 = struct('settings',st);
    cert1 = pielr_certify(PIE1,o1);
else
    fprintf('\n=== n = 1: shipped face, re-verified live ===\n');
    F = load(fullfile(here,'face_n1.mat'));      % struct 'face': V,S,Ns,part,...
    PIE1 = nb_rd2d(1);
    o1 = struct('settings',st,'face',F.face);
    cert1 = pielr_certify(PIE1,o1);
end
rows = {};
rows{end+1} = demo_row(1,cert1);
if ~cert1.ok
    error('demo_lowrank_2d: the n=1 certificate FAILED live verification; aborting.');
end
certs = struct();  certs.n1 = cert1;                                       % CC, 09/19/2026
if isempty(cert1.part)
    error('demo_lowrank_2d: no group partition at n=1; cannot tensor the face.');
end

% ---- the ladder: n >= 2 via the tensored face, no discovery -----------------
for n = nlist(nlist>=2)
    fprintf('\n=== n = %d: tensored n=1 face, certified live ===\n',n);
    faceN = pielr_tensor(cert1,n);
    PIEn = nb_rd2d(n);
    on = struct('settings',st,'face',faceN);
    certn = pielr_certify(PIEn,on);
    rows{end+1} = demo_row(n,certn); %#ok<SAGROW>
    print_table(rows);
    if ~certn.ok
        error('demo_lowrank_2d: the n=%d certificate FAILED live verification; aborting.',n);
    end
    certs.(sprintf('n%d',n)) = certn;                                   % CC, 09/19/2026
end

% ---- final table + reference lines ------------------------------------------
fprintf('\n==================== FINAL TABLE ====================\n');
print_table(rows);
fprintf(['\nReference (measured on this LPI family, same settings):\n' ...
    '  stock full solve of the SAME n=1 program: 1707.3 s (op rel 1.632e-07)\n' ...
    '  stock route at n=4: needs a 66.8 GiB dense normal matrix ->\n' ...
    '  memory-infeasible on a 63.7 GiB machine, while the face route above\n' ...
    '  certifies n=4 in minutes (dominated by LPI assembly, not solving).\n']);
if RUN_STOCK
    fprintf('\n=== stock full solve at n=1 (RUN_STOCK = true; ~28 min) ===\n'); %#ok<UNRCH> % config constant
    run_stock_n1(st);
end
fprintf('(certificates remain in your workspace: certs.n1, certs.n2, ...)\n'); % CC, 09/19/2026
fprintf('\nDEMO_LOWRANK_2D_DONE\n');

% ==================== script-local functions =================================
function row = demo_row(n,cert)
row = struct('n',n,'unk',cert.unknowns_full, ...
    'setup',cert.t_setup,'solve',cert.t_discover+cert.t_certify, ...
    'rel',cert.op_rel,'psd',cert.psd,'route',cert.route,'ok',cert.ok, ...
    'r',cert.r);
end

function print_table(rows)
fprintf('\n %2s | %-13s | %-8s | %-8s | %-11s | %-10s | %-3s | %s\n', ...
    'n','Gram unknowns','setup_s','solve_s','solve/setup','op_rel','psd','route');
for k = 1:numel(rows)
    r = rows{k};
    fprintf(' %2d | %13d | %8.1f | %8.1f | %11.4f | %10.3e | %3d | %s\n', ...
        r.n,r.unk,r.setup,r.solve,r.solve/max(r.setup,eps),r.rel,r.psd,r.route);
end
end

function run_stock_n1(st) %#ok<DEFNU> % called only when RUN_STOCK = true
% the stock reference: solve the full LPI with SeDuMi through lpisolve, as the
% 2-D executive would.  Kept behind a flag because it is ~28 min of pure solver
% time for a number the table above already quotes from measurement.
PIE = nb_rd2d(1);
t0 = tic;
prog = build_stab_2d_st2(PIE,st); % same program as the face route
sopts.solver = 'sedumi';  sopts.simplify = 0;
prog = sossolve(prog,sopts);
fprintf('stock full solve: %.1f s  (feasratio %.3g)\n', ...
    toc(t0),prog.solinfo.info.feasratio);
end

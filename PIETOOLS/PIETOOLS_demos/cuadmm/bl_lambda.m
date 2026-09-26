function L = bl_lambda(phys,setname,exec,lo,hi,rtol,tag,deadline)    % CC, 09/26/2026
% BL_LAMBDA  Locate a 1-D stability LPI's own boundary lambda_LPI by REBUILDING
% the program at each fraction of the PDE's analytic threshold, each program
% decided by bl_bisect in 'feas' mode with Mosek (verified F or Farkas I only).
%
%   L = bl_lambda('rd','heavy','PIE2PDEstability',0.95,1.0,1e-4,'lamhv')
%
% Why: the stability sentinels must sit just outside the LPI's boundary, not
% only outside the PDE's (lambda* = pi^2), because lambda_LPI < lambda* and a
% sentinel far past lambda_LPI tests only gross unsoundness of the verifier.
% hi is taken as infeasible by theorem (frac >= 1: the PDE is not stable), lo
% must be feasible (it is checked).  An undecided rebuild (Mosek neither F nor
% I) is treated as the infeasible side and recorded, so lam_I is then an upper
% estimate, not a certificate: L.lam_I_certified says which.
%
% Each rebuild goes through bl_run (struct-chunk form), the same dump path as
% the registry, so the program solved is exactly the one dumped.

if nargin < 8 || isempty(deadline), deadline = Inf; end
cuadmm_path;
DMP = fullfile(cuadmm_outdir(),'baseline','dumps');
L = struct('phys',phys,'set',setname,'exec',exec,'steps',[],'lam_F',NaN, ...
           'lam_I',hi,'lam_I_certified',false,'lam_I_id','');
a = lo;  b = hi;  first = true;
while (b - a) > rtol*b
    if posixtime(datetime('now')) + 60 > deadline, L.note = 'deadline'; break; end
    if first, f = lo; else, f = (a + b)/2; end
    id = sprintf('%s_f%.6f',tag,f);  id = strrep(id,'.','p');
    c  = struct('id',id,'cls','lambda','dim','1D','kind','feas', ...
                'builder','bl_b_stab1','args',{{phys,f,setname,exec}});
    bl_run(c);
    R  = bl_bisect(fullfile(DMP,[id '.mat']),struct('mode','feas','solver','mosek', ...
                   'tag','lam','deadline',deadline));
    v  = R.probes(end).verdict;
    L.steps(end+1,:) = [f, double(v=='F'), double(v=='I')];
    fprintf('LAMBDA %s frac=%.6f -> %s\n',tag,f,v);
    if first
        if v ~= 'F', L.note = sprintf('lo=%.4f is not certified feasible (%s)',lo,v); return; end
        a = f;  L.lam_F = f;  first = false;  continue
    end
    if v == 'F'
        a = f;  L.lam_F = f;
    else
        b = f;  L.lam_I = f;  L.lam_I_certified = (v == 'I');  L.lam_I_id = id;
    end
end
fprintf('LAMBDADONE %s  lambda_LPI in [%.6f, %.6f] (upper end certified: %d)\n', ...
        tag,L.lam_F,L.lam_I,L.lam_I_certified);
end

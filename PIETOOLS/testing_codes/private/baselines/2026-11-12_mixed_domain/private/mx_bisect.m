function B = mx_bisect(sys,par,base,lo,hi,nstep,tier)
% B = MX_BISECT(SYS,PAR,BASE,LO,HI,NSTEP,TIER) bisects the reach of
% PIETOOLS_stability_sop ('P') in one parameter of a system of MX_SYS.
% INPUT
% - sys:    char, the system name;
% - par:    char, the parameter bisected;
% - base:   'struct', the other parameters;
% - lo, hi: double, a certified and an uncertified value of PAR;
% - nstep:  integer, the number of bisection steps;
% - tier:   char, the lpisettings_sop tier.
% OUTPUT
% - B:      'struct' with fields lo, hi (the final bracket) and rec (the
%           MX_EXEC record of each step); one BISECT line per step.

% Initial coding MMP, 11/12/2026: bisection of N8 (README.md).
B = struct('lo',lo,'hi',hi,'rec',[]);
for k=1:nstep
    mid = (B.lo+B.hi)/2;
    p = base;  p.(par) = mid;
    R = mx_exec({sys,p,'stab','P',true},tier);
    B.rec = [B.rec; R];
    if R.certified, B.lo = mid; else, B.hi = mid; end
    fprintf('BISECT %s %s: step %d, %s = %.6g -> %s; bracket [%.6g, %.6g]\n',sys,tier,k,par,mid,R.status,B.lo,B.hi);
end
end

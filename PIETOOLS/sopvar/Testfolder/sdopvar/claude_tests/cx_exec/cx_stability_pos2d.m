function [prog,Pm] = cx_stability_pos2d(prog,X,side,dom,deg,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,PM] = CX_STABILITY_POS2D(PROG,X,SIDE,DOM,DEG,OPTS) is the container
% counterpart of the 2-D stock call  poslpivar_2d(prog,X.dim,DEG,OPTS)  for a
% state purely in L2[s1,s2]: a psd 'cdopvar' (poscopvar) on the SIDE
% ('out'|'in') space of the 1x1 container X.
%
% Degrees and options go through settings2possopvar, whose L2[s1,s2]
% translation is EXACT (full subset-cap arrays, exclude(8:16) -> include,
% sep(3:6) -> per-direction sep; settings2possopvar.m:29-45, 179-274) and
% whose vocabulary copquadvar shares for a space holding every registry
% variable (copquadvar.m:441-447). What cannot be reproduced errors here:
%   psatz other than 0/1 (copquadvar.m:347-350; settings2possopvar would
%   silently map 3..6 to 1 in block2d), anything the translator reports as
%   lossy (mixed sep patterns), and any R, L2[x] or L2[y] component (no
%   tested degree map; test_poscopvar_vs_poslpivar_2d.m:60-66).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sp,dm] = cx_space_list(X,side);
if numel(sp)~=1 || numel(sp{1})~=2
    error('cx_stability:scope2D',['Only a state purely in L2[s1,s2] is transcribed '...
          '(no tested R/L2[x]/L2[y] degree map for poscopvar).'])
end
pz = 0;     if isfield(opts,'psatz') && ~isempty(opts.psatz),   pz = opts.psatz;    end
if ~ismember(pz,[0 1])
    error('cx_stability:psatz','poslpivar_2d psatz=%g has no poscopvar form (0/1 only).',pz)
end
% One-term settings struct: LF and eq both set, so the translator reports no
% spurious 'missing field' loss for the half not used.
s = struct('is2D',1,'LF_deg',deg,'LF_opts',opts,'eq_deg',deg,'eq_opts',opts);
o = settings2possopvar(s);
if ~isempty(o.report.lossy)
    error('cx_stability:lossy','Not reproducible in poscopvar: %s',strjoin(o.report.lossy,' | '))
end
popt = struct('include',{{o.LF.opts.include}},'sep',o.LF.opts.sep,'psatz',pz);
[prog,Pm] = poscopvar(prog,dm,sp,dom,{o.LF.deg},popt);
end

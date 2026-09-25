function out = bl_fix(matfile,gam,outdir)
% bl_fix(matfile,gam,outdir) -- turn an objective-form SDP into the FEASIBILITY
% question at a fixed gamma, exactly, with no re-derivation of the LPI.
%
% WHY AT THE SDP LEVEL. The stock 1-D executives hard-declare the objective:
%   dpvar gam; lpidecvar(prog,gam); lpi_ineq(prog,gam); lpisetobj(prog,gam)
% (e.g. PIETOOLS_Hinf_gain.m:128-131), and only the 2-D executives accept a
% `gain` argument. Hand-mirroring eight executives to pin gamma would mean eight
% chances to build a subtly different program, and the mirror could not be
% checked against the original except by the numbers it produces. The objective
% vector is a single unit entry (verified in the dumps: C.txt reads "43 0 1"),
% so the same thing is achieved exactly by appending ONE row:
%       e_j' x = gam,      c := 0
% That is the bisection step by definition, and it is the identical SDP the
% executive would pose with gamma pinned.
%
% ORDER OF OPERATIONS MATTERS. sossolve normalises b, so b is rescaled AFTER the
% new row is appended -- appending to an already-normalised b would pin gamma to
% the wrong value by the old scale factor.
%
% Returns the new shape plus j and the un-normalised gam so a caller can undo the
% scaling when reading a solution back.

S = load(matfile);
c = S.c(:);
j = find(c);
if numel(j)~=1
    error('bl_fix:obj','expected exactly one objective entry, found %d',numel(j));
end
if abs(full(c(j))-1) > 1e-12
    error('bl_fix:coef','objective coefficient is %g, not 1',full(c(j)));
end
if j > S.K.f
    error('bl_fix:cone','objective index %d is not in the free block (K.f=%d)',j,S.K.f);
end

m    = numel(S.b);
nvar = size(S.At,1);
% b in the dump is already divided by bscl; recover the unnormalised b first so
% the appended row carries gam in the SAME units as the rest of the system.
Mt = load(strrep(matfile,'.mat','_meta.mat'));
b_un = full(S.b(:))*Mt.bscl;

row = sparse(j,1,1,nvar,1);
At2 = [S.At, row];
b2  = [b_un; gam];
bscl2 = norm(b2);  if bscl2==0 || ~isfinite(bscl2), bscl2 = 1; end

D.At = At2;
D.b  = b2/bscl2;
D.c  = sparse(nvar,1);              % feasibility: no objective at all
D.K  = S.K;  D.Ns = S.K.s;  D.Kf = S.K.f;

if nargin>=3 && ~isempty(outdir)
    if ~exist(fileparts(outdir),'dir'), mkdir(fileparts(outdir)); end
    save([outdir '.mat'],'-struct','D','-v7.3');
    RR = Mt.RR; bscl = bscl2; Sshape = Mt.S; obj_idx = j; gam_fixed = gam;
    save([outdir '_meta.mat'],'RR','bscl','Sshape','obj_idx','gam_fixed');
    dump2cuadmm([outdir '.mat'],outdir);
end

out = struct('m',m+1,'nvar',nvar,'obj_idx',j,'bscl',bscl2,'gam',gam, ...
             'Kf',S.K.f,'Ks',S.K.s,'At',{D.At},'b',{D.b},'c',{D.c},'K',D.K);
end

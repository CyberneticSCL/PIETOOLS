% fixednorm.m -- (a) assert the innerprod fix is the copy MATLAB actually loads,
% and (b) measure candidate FIXED operator-space normalisers.
%
% The problem: a gate normalised by ||Dop|| or ||Pop|| inherits their swing,
% because both are SOLUTION-dependent (measured: ||Dop|| varies 24x between two
% formulations of identical physics).  ||b|| is fixed but lives in coefficient
% space.  Wanted: something fixed, in OPERATOR space, same units as the residual.
%
% Candidate: the P = I negativity operator  N := T'A + A'T.  It is determined by
% the PLANT alone -- no decision variables, no settings degrees, no solution --
% and it is exactly the object the residual is a perturbation of.
cuadmm_path;
fprintf('FN which(innerprod) = %s\n', which('innerprod'));
fprintf('FN n_resolutions = %d\n', numel(unique(which('innerprod','-all'))));
fprintf('FN first line = %s\n', strtrim(fileread_line(which('innerprod'))));

pvar s t
for lam = [0 2 5]
    x = pde_var('state',1,s,[0,1]);
    PIE = initialize(convert([diff(x,t,1)==diff(x,s,2)+lam*x;
          subs(x,s,0)==0; subs(x,s,1)==0]));
    Top = PIE.T;  Aop = PIE.A;
    opvar Iop; Iop.I=PIE.dom; Iop.var1=PIE.vars(1,1); Iop.var2=PIE.vars(1,2);
    Iop.dim = Top.dim;
    if Top.dim(1,1)>0, Iop.P = eye(Top.dim(1,1)); end
    if Top.dim(2,1)>0, Iop.R.R0 = eye(Top.dim(2,1)); end
    N  = clean_opvar(Aop'*(Iop*Top) + (Iop*Top)'*Aop, 1e-12);   % P = I negativity
    TT = Top'*Top;
    nN  = norm(cuadmm_private('opnorm_pi',N ,20));
    nTT = norm(cuadmm_private('opnorm_pi',TT,20));
    nT  = norm(cuadmm_private('opnorm_pi',Top,20));
    nA  = norm(cuadmm_private('opnorm_pi',Aop,20));
    fprintf(['FN lam=%-3g  ||T''A+A''T||=%.6e  ||T''T||=%.6e  ||T||=%.6e  ' ...
             '||A||=%.6e  ||T||*||A||=%.6e\n'], lam, nN, nTT, nT, nA, nT*nA);
end
fprintf('FNDONE\n');

function l = fileread_line(f)
fid = fopen(f,'r'); l = fgetl(fid); fclose(fid);
end

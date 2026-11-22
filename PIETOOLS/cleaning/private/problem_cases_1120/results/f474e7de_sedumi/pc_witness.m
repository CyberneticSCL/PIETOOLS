function pc_witness()
% PC_WITNESS() measures whether the positive "numerical spectrum max Re" that
% PIETOOLS_stability_sop prints for the 1-D position forms of DIAG44F and
% DIAG69 is a property of the system or of the discretization: the largest
% real part of the finite eigenvalues of the Legendre-Galerkin pencil (A,T)
% (pie_witness_sop, mode 'rate') at N = 12, 24 (the default), 48 and 96.
% Systems: 'dw1' (u_tt = u_ss, u(0) = 0, u_s(1) = -u_t(1); no eigenvalues),
% '69pos' (DIAG69), '44min' (DIAG44E), and the energy form of PDE 68 as a
% control. results_WITNESS.csv in this folder. A diagnostic.
% MMP, 11/20/2026.
here = fileparts(mfilename('fullpath'));    cd(here);
of = fullfile(here,'results_WITNESS.csv');
fid = fopen(of,'w');  fprintf(fid,'system,N,ok,maxre,note\n');  fclose(fid);
for v = {'dw1','69pos','44min','68'}
    PIE = buildw(v{1});
    for N = [12 24 48 96]
        try
            W = pie_witness_sop(PIE,'rate',struct('N_cheb',N));
            mr = NaN;   if W.ok,  mr = W.maxre;  end
            nt = '';    if isfield(W,'note'),  nt = strrep(W.note,',',';');  end
            ok = W.ok;
        catch ME
            ok = false;  mr = NaN;  nt = strrep(ME.message,',',';');
        end
        fid = fopen(of,'a');  fprintf(fid,'%s,%d,%d,%.6g,"%s"\n',v{1},N,ok,mr,nt);  fclose(fid);
        fprintf('##### WITNESS %s N %d: ok %d maxre %.6g %s\n',v{1},N,ok,mr,nt);
    end
end
fprintf('## WITNESS DONE %s\n',char(datetime('now')));
end

function PIE = buildw(v)
pvar s t
clear stateNameGenerator
switch v
    case 'dw1'
        u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);
        PDE_t = [diff(u,'t')==v1;  diff(v1,'t')==diff(u,s,2);
                 subs(u,s,0)==0;  subs(v1,s,0)==0;  subs(diff(u,s),s,1)==-subs(v1,s,1)];
    case '69pos'
        u = pde_var(s,[0,1],2);  vv = pde_var(s,[0,1],1);  w = pde_var(s,[0,1],2);
        PDE_t = [diff(u,'t')==vv;  diff(vv,'t')==diff(u,s,2);  diff(w,'t')==diff(w,s,2);
                 subs(u,s,1)==0;  subs(diff(u,s),s,0)==-subs(diff(w,s),s,0);
                 subs(vv,s,1)==0;  subs(w,s,1)==0;  subs(w,s,0)==subs(vv,s,0)];
    case '44min'
        k = 1;  tau = 1;  mu = 0.4;
        u = pde_var(s,[0,1],2);  v1 = pde_var(s,[0,1],1);  z = pde_var(s,[0,1],1);
        PDE_t = [diff(u,'t')==v1;  diff(v1,'t')==diff(u,s,2);  diff(z,'t')==-(1/tau)*diff(z,s);
                 subs(u,s,0)==0;  subs(v1,s,0)==0;
                 subs(diff(u,s),s,1)==-k*(1-mu)*subs(v1,s,1)-k*mu*subs(z,s,1);
                 subs(z,s,0)==subs(v1,s,1)];
    case '68'
        PDE_t = PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty(0,{});
end
PIE = convert(PDE_t,'pie');
end

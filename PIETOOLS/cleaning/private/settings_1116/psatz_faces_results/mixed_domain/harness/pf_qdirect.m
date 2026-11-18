function R = pf_qdirect(sysname,plist,qdegs,dws,tag,eppos,epneg,topt)
% R = PF_QDIRECT(SYS,PLIST,QDEGS,DWS,TAG,EPPOS,EPNEG) the PIE2PDE condition
% in its direct form (Cor. 7.1 of the hierarchy paper): a free PI operator
% Q with the shape of T (LPIVAR_CDOPVAR at the degrees QDEGS, a cell of
% degree specs), the symmetry T'*Q = Q'*T as an equality, the coercivity
% T'*Q >= EPPOS T'*T and the decay -(A'*Q + Q'*A) - 2 EPNEG T'*Q >= 0 as
% inequalities certified by LPI_INEQ_SOP, each sized by the reader from
% its own support (light tier: dD 1, dw = DW, terms 'auto'). No storage
% operator. Results in results/qdirect_<TAG>.mat. SCRATCH, MMP, 11/16/2026.
if nargin<6 || isempty(eppos),  eppos = 1e-4;   end
if nargin<7 || isempty(epneg),  epneg = 0;      end
if nargin<8 || isempty(topt),   topt = struct(); end    % term options merged into both inequalities (psatz_rows etc.)
S = pf_setup(false);
out = sprintf('%s/results/qdirect_%s.mat',S,tag);
R = struct('sys',{},'par',{},'qdeg',{},'terms',{},'dw',{},'status',{},'ndv',{},'t',{},'numerr',{},'feasratio',{},'relres',{},'eta',{},'err',{});
for ip = 1:numel(plist)
    p = plist{ip};
    Sys = mixed_domain_systems(sysname,p);
    PIE = exec_pie_sop(initialize(convert(Sys.PDE,{'pie','silent'})));
    Tm = PIE.T;     Am = PIE.A;
    for iq = 1:numel(qdegs)
        qd = qdegs{iq};
        for dw = dws
            r = struct('sys',sysname,'par',p,'qdeg',jsonencode(qd),'terms',jsonencode(topt),'dw',dw,'status','','ndv',NaN,'t',NaN, ...
                       'numerr',NaN,'feasratio',NaN,'relres',NaN,'eta',NaN,'err','');
            t0 = tic;
            try
                st = lpisettings_sop('light');  st.sos_opts.solver = 'mosek';
                prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);
                [prog,Qm] = lpivar_cdopvar(prog,Tm,qd);
                prog = lpi_eq_sop(prog,Tm'*Qm - Qm'*Tm);
                Pq = 0.5*(Tm'*Qm + Qm'*Tm);
                oi = topt;  oi.dw = dw;
                [prog,~,~] = lpi_ineq_sop(prog,Pq - eppos*(Tm'*Tm),oi);
                Dm = Am'*Qm + Qm'*Am + 2*epneg*Pq;
                [prog,~,~] = lpi_ineq_sop(prog,-Dm,oi);
                r.ndv = size(prog.decvartable,1);
                prog = lpisolve(prog,st.sos_opts);
                cl = lpi_classify_sop(prog,[]);
                r.t = toc(t0);  r.status = cl.status;   r.numerr = cl.numerr;
                r.feasratio = cl.feasratio;     r.relres = cl.relres;  r.eta = cl.eta;
            catch ME
                r.err = strrep(ME.message,newline,' ');     r.t = toc(t0);
            end
            R(end+1) = r;                                                   %#ok<AGROW>
            fprintf('@@ QD %s par=%s qdeg=%s terms=%s dw=%d %s ndv=%d t=%.0f numerr=%g feas=%.3g relres=%.1e eta=%.1e %s\n', ...
                    sysname,jsonencode(p),r.qdeg,r.terms,dw,r.status,r.ndv,r.t,r.numerr,r.feasratio,r.relres,r.eta,r.err);
            save(out,'R');
        end
    end
end
fprintf('## QD DONE %s\n',tag);
end

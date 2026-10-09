function [prog,P,R,omega,info] = PIETOOLS_well_posedness_sop(PIE,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,R,OMEGA,INFO] = PIETOOLS_WELL_POSEDNESS_SOP(PIE,SETTINGS) the
% container version of PIETOOLS_well_posedness: existence of P >= eppos I,
% R >= I with
%     T'PA + A'PT <= 2 omega T'PT      (dissipativity)
%     (T - A) R (T - A)' >= eppos I    (surjectivity)
% omega = -settings.epneg (positive allowed), which gives ||e^{tA}|| <= M
% e^{omega t}. The storage P, R are the stock ones (lf degrees, 'out' and
% 'in' of T), the two inequalities LPI_INEQ_SOP with 'like' P and 'like' R
% (weight raised by settings.sop.dw_P), or the stock slack when
% settings.sop.slack is 'stock'; rebuilt while the solver returns UNKNOWN or
% INFEASIBLE (LPI_SOLVE_LOOP_SOP). Without SETTINGS the stock defaults:
% heavy, eppos = eppos2 = 1e-2, epneg = 0. 1-D only, as the stock.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2 || isempty(settings)
    settings = lpisettings('heavy');
    settings.sos_opts.simplify = 1;
    settings.eppos = 1e-2;  settings.eppos2 = 1e-2;   settings.epneg = 0;
end
label = 'well_posedness_sop';
st = exec_sop_settings(settings);
if ~isa(PIE,'pie_struct'),  error('PIETOOLS_well_posedness_sop:pie','The PIE should be a ''pie_struct''.');  end
PIE = initialize(PIE);
if PIE.dim==2,  error('PIETOOLS_well_posedness_sop:dim','Well-posedness analysis of 2D PDEs is currently not supported.');   end
tl = exec_tools_sop(PIE,st);
eppos = tl.eppos;   omega = -tl.epneg;
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
P = [];     R = [];
try
    P = tl.toop(getsol_lpivar_sop(prog,aux.Pm));
    R = tl.toop(getsol_lpivar_sop(prog,aux.Rm));
catch ME
    if st.sop.verbose,  fprintf('  [%s] solution extraction: %s\n',label,ME.message);  end
end
info = exec_info_sop(prog,aux,st,label);
info.omega = omega;     info.eppos = eppos;
info.D1 = aux.Dm1;      info.D2 = aux.Dm2;
info.witness = [];      info.maxre = NaN;
if st.sop.witness
    Wt = pie_witness_sop(PIE,'rate',st.sop);
    info.witness = Wt;
    if Wt.ok,   info.maxre = Wt.maxre;  end
end
if st.sop.verbose
    fprintf('  [%s] %s: certificate at omega %g: %d; numerical spectrum max Re %.6g\n', ...
            label,cl.status,omega,info.certified,info.maxre);
end

    function [prog,aux] = assemble(dw)
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end     % ask lpi_eq_cdopvar to record its rows
        [prog,Pm] = tl.storage(prog,Tm,'out');
        [prog,Rm] = tl.storage(prog,Tm,'in');
        I1 = tl.eye(1,eppos);
        Pm = Pm + I1;                               % P >= eppos I
        Rm = Rm + tl.eye(2,[1;1]);                  % R >= I
        Dm1 = Tm'*Pm*Am + Am'*Pm*Tm;
        if omega~=0,    Dm1 = Dm1 - 2*omega*(Tm'*Pm*Tm);    end
        TA = Tm - Am;
        Dm2 = TA*Rm*TA' - I1;
        if strcmp(st.sop.slack,'stock')
            [prog,N1] = tl.slack(prog,Dm1);
            [prog,tag] = lpi_eq_sop(prog,Dm1+N1,'symmetric');
            [prog,N2] = tl.slack(prog,Dm2);
            prog = lpi_eq_sop(prog,Dm2-N2,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            o1 = struct('like',Pm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin);
            [prog,~,li] = lpi_ineq_sop(prog,-Dm1,o1);
            o2 = struct('like',Rm,'dw',st.sop.dw_P+dw,'margin',st.sop.margin);
            [prog,~,li2] = lpi_ineq_sop(prog,Dm2,o2);
            li.degrees2 = li2.degrees;
        end
        aux = struct('obj',[],'Pm',Pm,'Rm',Rm,'Dm1',Dm1,'Dm2',Dm2,'li',li);
    end
end

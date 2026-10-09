function [prog,Pout,info] = stability_build_sop(PIE,st,form,label)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P,INFO] = STABILITY_BUILD_SOP(PIE,ST,FORM,LABEL) the stability
% executives on the container path, 1-D and 2-D, one LPI per FORM:
%
%   'P'   PIETOOLS_PDEstability, _stability_2D:   P >= 0 (+ eppos I),
%         T'PA + A'PT + epneg T'PT <= 0 (2-D: 2 epneg T'PT, as stock);
%   'Pd'  PIETOOLS_PDEstability_dual, _stability_dual_2D:
%         T P A' + A P T' + epneg T P T' <= 0 (2-D: 2 epneg);
%   'Q'   PIETOOLS_PIE2PDEstability (1-D):   P >= 0 (+ eppos2 T'T), T'Q = P,
%         A'Q + Q'A + epneg P <= 0;
%   'Qd'  PIETOOLS_PIE2PDEstability_dual (1-D): P + eppos2 T T', TQ = P,
%         AQ + Q'A' + epneg P <= 0.
%
% The storage is the stock one (EXEC_TOOLS_SOP storage, 'out' or 'in' as the
% stock; eppos the per-space weights of the stock: eppos, eppos2 in 1-D, the
% four-entry eppos of settings_2d in 2-D), the negativity constraint is
% LPI_INEQ_SOP with 'like' the storage (weight raised by ST.sop.dw_P in the
% coercive forms, dw_Q in the Q forms), or the stock slack when ST.sop.slack
% is 'stock'; the program is rebuilt while the solver returns UNKNOWN or
% INFEASIBLE (LPI_SOLVE_LOOP_SOP). A feasibility LPI: INFO.certified is true
% when the solver found the certificate at the given epneg. In 1-D
% INFO.witness is the numerical spectrum of (A,T) (PIE_WITNESS_SOP 'rate'):
% INFO.maxre, the largest real part, to set beside the certified decay.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<4,    label = ['stability ' form];    end
st = exec_sop_settings(st);
if ~isa(PIE,'pie_struct'),  error('stability_build_sop:pie','The PIE should be a ''pie_struct''.');     end
PIE = initialize(PIE);
tl = exec_tools_sop(PIE,st);    nv = tl.nv;
if ~any(strcmp(form,{'P','Pd','Q','Qd'}))
    error('stability_build_sop:form','FORM should be ''P'', ''Pd'', ''Q'' or ''Qd''.');
end
if nv==2 && any(strcmp(form,{'Q','Qd'}))
    error('stability_build_sop:Q2D','The Q forms have no stock 2-D counterpart; PIETOOLS_stability_2D_sop and _dual_2D_sop are the coercive forms.');
end
coercive = any(strcmp(form,{'P','Pd'}));
dwbase = st.sop.dw_Q;   if coercive,    dwbase = st.sop.dw_P;   end
eppos = tl.eppos;   epneg = tl.epneg;
epfac = 1;  if nv==2,   epfac = 2;  end           % the stock 2-D executives write 2 epneg
Tm = op2copvar_sop(PIE.T);  Am = op2copvar_sop(PIE.A);
[prog,aux] = lpi_solve_loop_sop(@assemble,st,label);
cl = aux.classify;
Psol = [];  Pout = [];
try
    Psol = getsol_lpivar_sop(prog,aux.Pm);
    Pout = tl.toop(Psol);
catch
    Pout = Psol;
end
info = exec_info_sop(prog,aux,st,label);
info.form = form;   info.epneg = epneg;     info.eppos = eppos;     info.nv = nv;
info.storage = Psol;    info.D = aux.Dm;    info.N = aux.Nm;
info.witness = [];  info.maxre = NaN;
if st.sop.witness && nv==1
    Wt = pie_witness_sop(PIE,'rate',st.sop);
    info.witness = Wt;
    if Wt.ok,   info.maxre = Wt.maxre;  end
end
if st.sop.verbose
    fprintf('  [%s] %s: certificate at epneg %g: %d; numerical spectrum max Re %.6g\n', ...
            label,cl.status,epneg,info.certified,info.maxre);
end

    function [prog,aux] = assemble(dw)
        prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
        if st.sop.keepdual,     prog.sopeq = {};    end     % ask lpi_eq_cdopvar to record its rows
        switch form
            case 'P'
                [prog,Pm] = tl.storage(prog,Tm,'out');
                if ~all(eppos==0),  Pm = Pm + tl.eye(1,eppos);  end
                Dm = Tm'*Pm*Am + Am'*Pm*Tm;
                if epneg~=0,    Dm = Dm + epfac*epneg*(Tm'*Pm*Tm);  end
            case 'Pd'
                [prog,Pm] = tl.storage(prog,Tm,'in');
                if ~all(eppos==0),  Pm = Pm + tl.eye(2,eppos);  end
                Dm = Tm*Pm*Am' + Am*Pm*Tm';
                if epneg~=0,    Dm = Dm + epfac*epneg*(Tm*Pm*Tm');  end
            case 'Q'
                [prog,Pm] = tl.storage(prog,Tm,'in');
                Pm = Pm + eppos(2)*(Tm'*Tm);
                [spO,dmO] = copvar_space_list(Tm,'out');    [spI,dmI] = copvar_space_list(Tm,'in');
                [prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pm));
                prog = lpi_eq_sop(prog,Tm'*Qm - Pm);
                Dm = Am'*Qm + Qm'*Am;
                if epneg~=0,    Dm = Dm + epneg*Pm;     end
            case 'Qd'
                [prog,Pm] = tl.storage(prog,Tm,'out');
                Pm = Pm + eppos(2)*(Tm*Tm');
                [spO,dmO] = copvar_space_list(Tm,'in');     [spI,dmI] = copvar_space_list(Tm,'out');
                [prog,Qm] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pm));
                prog = lpi_eq_sop(prog,Tm*Qm - Pm);
                Dm = Am*Qm + Qm'*Am';
                if epneg~=0,    Dm = Dm + epneg*Pm;     end
        end
        if strcmp(st.sop.slack,'stock')
            [prog,Nm] = tl.slack(prog,Dm);
            [prog,tag] = lpi_eq_sop(prog,Dm+Nm,'symmetric');
            li = struct('tag',tag,'degrees',[],'terms',[]);
        else
            o = struct('like',Pm,'dw',dwbase+dw,'margin',st.sop.margin);
            [prog,Nm,li] = lpi_ineq_sop(prog,-Dm,o);
        end
        aux = struct('obj',[],'Pm',Pm,'Dm',Dm,'Nm',Nm,'li',li,'form',form);
    end
end

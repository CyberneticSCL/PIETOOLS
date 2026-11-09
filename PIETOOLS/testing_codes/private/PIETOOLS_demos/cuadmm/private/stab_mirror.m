function out = stab_mirror(PIE,st)
% stab_mirror -- PIETOOLS_stability_sop in form 'Q' (the PIE2PDE form), built as
% that executive builds its first pass but NOT solved, RETURNING the operators.
% bl_scale and bl_big need the SDP's shape before deciding to solve (Mosek's
% dense Schur complement), and the executive solves internally.
%
% FAITHFULNESS IS CHECKED, NOT ASSUMED: sdpshape(out.prog) must equal sdpshape
% of the executive's solved program (m, Kf, Ns, nnz) on the same PIE and
% settings.
%
% CC, 11/07/2026: ported to the container executive. Deleted: the mirror of the
%   stock PIETOOLS_PIE2PDEstability as of 2026-09-23 (poslpivar for P,
%   get_lpivar_degs and lpivar for Q, lpi_eq for T'*Q = P, lpi_ineq or a
%   poslpivar slack with lpi_eq for the negativity), which stops on Top.dim
%   since convert returns copvar operators (3df33dcc, 10/29/2026); the opvar
%   removal retires that executive. Now the 'Q' branch of the executive's
%   assemble (MMP, 11/04/2026) at raise 0, the program it solves first, and
%   only with bl_settings' max_raise 0. 1-D only: the executive poses 'Q' as
%   'P' in two or more variables. Operators are cdopvar / copvar; out.Deop is
%   the slack N, Dop + Deop == 0 as before.
% MMP, 11/09/2026: the decay term 2*epneg*P, as the executive's 'Q' form
%   poses it since 11/09/2026 (one decay rate in every dimension).

if ~isa(PIE,'pie_struct'), error('stab_mirror:type','PIE must be a pie_struct'); end
% CC, 11/07/2026 (start): the executive's setup, then its 'Q' assemble at dw 0
PIE = exec_pie_sop(PIE);                    % every operator a copvar
st  = exec_sop_settings(st);                % container options completed
tl  = exec_tools_sop(PIE,st);   st.sop = tl.sop;
if tl.nv~=1, error('stab_mirror:nv','1-D only: %d variables pose ''Q'' as ''P''',tl.nv); end
dwbase = st.sop.dw_Q;                       % weight raise of the 'like' slack
if ~isempty(st.sop.dw_feas),    dwbase = st.sop.dw_feas(min(tl.nv,end));    end
eppos = tl.eppos;   epneg = tl.epneg;
if ~isempty(st.sop.eppos_stab), eppos(:) = st.sop.eppos_stab;   end
Top = PIE.T;    Aop = PIE.A;

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
if st.sop.keepdual, prog.sopeq = {};    end
% P >= eppos2*T'*T, Q free with T'*Q = P
[prog,Pop] = tl.storage(prog,Top,'in');
Pop = Pop + eppos(2)*(Top'*Top);
[spO,dmO] = copvar_space_list(Top,'out');   [spI,dmI] = copvar_space_list(Top,'in');
[prog,Qop] = lpivar_cdopvar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs_sop(Pop));
prog = lpi_eq_sop(prog,Top'*Qop - Pop);
% A'*Q + Q'*A + 2*epneg*P <= 0, by the slack the settings select
Dop = Aop'*Qop + Qop'*Aop;
% if epneg~=0,    Dop = Dop + epneg*Pop;  end                               % MMP, 11/09/2026 (was)
if epneg~=0,    Dop = Dop + 2*epneg*Pop;    end                             % MMP, 11/09/2026
if strcmp(st.sop.slack,'stock')
    [prog,Deop] = tl.slack(prog,Dop);
    prog = lpi_eq_sop(prog,Dop+Deop,'symmetric');
else
    o = struct('like',Pop,'dw',dwbase,'margin',st.sop.margin,'cap',tl.capfor(Dop));
    [prog,Deop] = lpi_ineq_sop(prog,-Dop,o);
end
% CC, 11/07/2026 (end)

out = struct('prog',prog,'Pop',Pop,'Qop',Qop,'Dop',Dop,'Deop',Deop,'Top',Top,'Aop',Aop);
end

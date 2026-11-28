function out = stab_mirror(PIE,st)
% stab_mirror -- PIETOOLS_stability in form 'Q' (the PIE2PDE form), built as
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
% CC, 11/17/2026: the executive's 'Q' branch is the direct form since
%   11/17/2026 (MMP; settings report 11/16/2026, Sec. 8.3): Q free with the
%   shape of T at lpi.var.deg, T'*Q - Q'*T = 0, P = sym(T'*Q) >= eppos2*T'*T
%   and the negativity, each by a slack sized from its own support (no
%   'like'), weight base 0; lpiprogram takes st (executive, 11/12/2026).
%   Mirrored. Before: light n=1 m 116 [10 14 7] here, 105 [16 8 18 9] in the
%   executive; heavy n=1, 2, 4 likewise. Old lines kept, marked (was).
% CC, 11/18/2026: Q at sop.deg_Q (lpi.var.deg when empty), the entry of the
%   PIE's dimension, as the executive's 'Q' branch since 11/18/2026 (MMP);
%   light's 1-D deg_Q is [4 2 2], lpi.var.deg [2 2 2]. First pass: no raise
%   of deg_Q or of the slack weight.

if ~isa(PIE,'pie_struct'), error('stab_mirror:type','PIE must be a pie_struct'); end
% CC, 11/07/2026 (start): the executive's setup, then its 'Q' assemble at dw 0
PIE = exec_pie(PIE);                        % every operator an opvar
st  = exec_sop_settings(st);                % container options completed
tl  = exec_tools(PIE,st);   st.sop = tl.sop;
% if tl.nv~=1, error('stab_mirror:nv','1-D only: %d variables pose ''Q'' as ''P''',tl.nv); end % CC, 11/17/2026 (was)
if tl.nv~=1, error('stab_mirror:nv','1-D only (checked in 1-D): %d variables',tl.nv); end % CC, 11/17/2026
% dwbase = st.sop.dw_Q;                       % weight raise of the 'like' slack % CC, 11/17/2026 (was)
% if ~isempty(st.sop.dw_feas),    dwbase = st.sop.dw_feas(min(tl.nv,end));    end % CC, 11/17/2026 (was)
dwbase = 0;     % 'Q': own-support slacks, base 0 as in the executive       % CC, 11/17/2026
eppos = tl.eppos;   epneg = tl.epneg;
if ~isempty(st.sop.eppos_stab), eppos(:) = st.sop.eppos_stab;   end
Top = PIE.T;    Aop = PIE.A;

% prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);                   % CC, 11/17/2026 (was)
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom,st);                  % CC, 11/17/2026
if st.sop.keepdual, prog.sopeq = {};    end
% CC, 11/17/2026 (start): the direct form. Was: storage P on the input
% spaces of T, Q at the degrees of P, T'*Q = P (lines below).
% % P >= eppos2*T'*T, Q free with T'*Q = P
% [prog,Pop] = tl.storage(prog,Top,'in');
% Pop = Pop + eppos(2)*(Top'*Top);
% [spO,dmO] = opvar_space_list(Top,'out');   [spI,dmI] = opvar_space_list(Top,'in');
% [prog,Qop] = lpivar(prog,struct('out',dmO,'in',dmI),struct('out',{spO},'in',{spI}),tl.dom,get_lpivar_degs(Pop));
% prog = lpi_eq(prog,Top'*Qop - Pop);
% Q free with the shape of T; T'*Q self-adjoint, P its symmetric part
% [prog,Qop] = lpivar_cdopvar(prog,Top,st.lpi.var.deg{min(tl.nv,end)});     % CC, 11/18/2026 (was); CC, 11/17/2026
qdeg = st.sop.deg_Q;    if isempty(qdeg),   qdeg = st.lpi.var.deg;  end     % CC, 11/18/2026
if ~iscell(qdeg),   qdeg = {qdeg};  end                                     % CC, 11/18/2026
[prog,Qop] = lpivar(prog,Top,qdeg{min(tl.nv,end)});                         % CC, 11/18/2026
prog = lpi_eq(prog,Top'*Qop - Qop'*Top);
Pop = 0.5*(Top'*Qop + Qop'*Top);
% P >= eppos2*T'*T, by a slack sized from its own support
o = struct('dw',dwbase,'margin',st.sop.margin,'cap',tl.capfor(Pop));
prog = lpi_ineq(prog,Pop - min(eppos(min(2,end):end))*(Top'*Top),o);
% CC, 11/17/2026 (end)
% A'*Q + Q'*A + 2*epneg*P <= 0, by the slack the settings select
Dop = Aop'*Qop + Qop'*Aop;
% if epneg~=0,    Dop = Dop + epneg*Pop;  end                               % MMP, 11/09/2026 (was)
if epneg~=0,    Dop = Dop + 2*epneg*Pop;    end                             % MMP, 11/09/2026
if strcmp(st.sop.slack,'stock')
    [prog,Deop] = tl.slack(prog,Dop);
    prog = lpi_eq(prog,Dop+Deop,'symmetric');
else
%   o = struct('like',Pop,'dw',dwbase,'margin',st.sop.margin,'cap',tl.capfor(Dop)); % CC, 11/17/2026 (was)
    o = struct('dw',dwbase,'margin',st.sop.margin,'cap',tl.capfor(Dop));    % CC, 11/17/2026
    [prog,Deop] = lpi_ineq(prog,-Dop,o);
end
% CC, 11/07/2026 (end)

out = struct('prog',prog,'Pop',Pop,'Qop',Qop,'Dop',Dop,'Deop',Deop,'Top',Top,'Aop',Aop);
end

function out = stab_mirror(PIE,st)
% stab_mirror -- PIETOOLS_PIE2PDEstability, reproduced line for line up to the
% solve, but RETURNING the operators so a residual can be measured at the
% operator level.  The stock executive returns only (prog,P), and the LPI's
% residual operator -- the thing that must be zero -- is internal to it.
%
% FAITHFULNESS IS CHECKED, NOT ASSUMED: t4 compares sdpshape(out.prog) against
% sdpshape of the stock executive's captured program.  If m, Kf, Ns and nnz all
% agree, the mirror poses the same SDP.
%
% Mirrors executives/PIETOOLS_PIE2PDEstability.m as of 2026-09-23.

if ~isa(PIE,'pie_struct'), error('stab_mirror:type','PIE must be a pie_struct'); end
PIE = initialize(PIE);
Top = PIE.T;   Aop = PIE.A;

dd1=st.dd1; dd12=st.dd12; options1=st.options1; options12=st.options12;
override1=st.override1; eppos=st.eppos; epneg=st.epneg; eppos2=st.eppos2;  %#ok<NASGU>

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);

% STEP 1: positive Lyapunov operator
[prog,P1op] = poslpivar(prog,Top.dim,dd1,options1);
if override1~=1
    [prog,P2op] = poslpivar(prog,Top.dim,dd12,options12);
    Pop = P1op+P2op;
else
    Pop = P1op;
end
Pop = Pop + eppos2*Top'*Top;

% Indefinite Qop with Top'*Qop = Pop
Qdeg = get_lpivar_degs(Pop,Top);
[prog,Qop] = lpivar(prog,Top.dim,Qdeg);
prog = lpi_eq(prog,Top'*Qop-Pop);

% STEP 2: the negativity operator
Dop = Aop'*Qop+Qop'*Aop + epneg*Pop;

% STEP 3: enforce negativity through an equality with a positive slack
if st.sosineq_on
    prog = lpi_ineq(prog,-Dop,st.opts);
    Deop = [];
else
    [prog,De1op] = poslpivar(prog,Dop.dim,st.dd2,st.options2);
    if st.override2~=1
        [prog,De2op] = poslpivar(prog,Dop.dim,st.dd3,st.options3);
        Deop = De1op+De2op;
    else
        Deop = De1op;
    end
    prog = lpi_eq(prog,Dop+Deop,'symmetric');
end

out = struct('prog',prog,'Pop',Pop,'Qop',Qop,'Dop',Dop,'Deop',Deop,'Top',Top,'Aop',Aop);
end

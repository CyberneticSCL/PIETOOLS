function [prog,H] = build_stab_st(PIE,st)
% PROVENANCE.  scratchpad/reach1d/build_stab_st.m verbatim.
%
% build_stab_st -- build_stability2 taking a settings STRUCT instead of a name.
%
% build_stability2 hardcodes st = lpisettings(setting), so it only accepts the
% stock names ('heavy','veryheavy',...).  The degree axis needs the numeric
% four-knob specs that lset() produces ([4 4 4 2] etc.), which is the whole
% point of that axis: it grows the Gram blocks WITHOUT adding states, and the
% banked minimum-rank study found minimum rank to be FLAT IN DEGREE (7,7,6 while
% the block grew 30 -> 33 -> 108).  Body is otherwise identical to
% build_stability2 so the LPI being solved is the same one.
% HARNESS COST, NOT SDP COST.  Requesting poslpivar's THIRD output (Qmat) makes
% it assemble the Gram as a dpvar; at spec [8 8 8 4] with n=2 that assembly asks
% for a 76.2 GB array and dies, and it dies sooner at [10 10 10 5].  Nothing here
% ever read H.Q -- it was a debugging convenience -- so all four calls now take
% two outputs and the high-degree rungs become reachable.  H.Q is left {}.
Top = PIE.T;    Aop = PIE.A;
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,P1op] = poslpivar(prog,Top.dim,st.dd1,st.options1);
Pop = P1op;
if st.override1~=1
    [prog,P2op] = poslpivar(prog,Top.dim,st.dd12,st.options12);
    Pop = P1op+P2op;
end
Imat = blkdiag(st.eppos*eye(Pop.dim(1,:)),st.eppos2*eye(Pop.dim(2,:)));
Pop = Pop + mat2opvar(Imat,Pop.dim(:,2),PIE.vars,PIE.dom);
Dop = Top'*Pop*Aop + Aop'*Pop*Top + st.epneg*Top'*Pop*Top;
[prog,De1op] = poslpivar(prog,Dop.dim,st.dd2,st.options2);
Deop = De1op;
if st.override2~=1
    [prog,De2op] = poslpivar(prog,Dop.dim,st.dd3,st.options3);
    Deop = De1op+De2op;
end
prog = lpi_eq(prog,Dop+Deop,'symmetric');
H.Top=Top; H.Aop=Aop; H.Pop=Pop; H.Dop=Dop; H.Deop=Deop;
H.Q={}; H.st=st; H.PIE=PIE;
end

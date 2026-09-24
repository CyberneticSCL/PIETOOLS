function [prog,H] = build_stab_1d(PIE,st,dual)                              % CC, 09/23/2026
% BUILD_STAB_1D  1-D direct-form stability LPI, returning OPERATOR HANDLES.
%
% The 1-D counterpart of build_stab_2d_st2, and the same LPI as
% executives/PIETOOLS_PDEstability.m (dual = false) and
% executives/PIETOOLS_PDEstability_dual.m (dual = true), minus the solve.
% Those two executives differ in exactly one line -- the negativity operator
% is T*P*A' + A*P*T' rather than T'*P*A + A'*P*T -- so one flag covers both.
%
% WHY THIS EXISTS.  The executive returns prog alone, and the operator handles
% are what the operator-level gate needs: with norm(b) ~ 5e-6 on these
% programs an equality residual cannot establish feasibility, so verification
% has to push the candidate back through Pop and Deop.  Same reasoning as
% build_stab_2d_st2's header, one dimension down.
%
% RELATION TO tests_1d/build_stab_st.m.  That file is the legacy 1-D fixture
% and builds the identical program; it is left alone so it can serve as the
% behaviour baseline.  T-neutral in the regression suite asserts that this
% routine reproduces its Atf, bf and block layout exactly.
%
% INPUT   PIE  pie_struct (or any struct with T, A, vars, dom), PIE.dim == 1
%         st   1-D settings struct, as lpisettings('light'|'heavy'|...)
%              returns: dd1, dd12, options1, options12, override1, eppos,
%              eppos2, epneg, override2, options2, options3, dd2, dd3
%         dual optional, default false
% OUTPUT  prog unsolved LPI program
%         H    .Top .Aop .Pop .Dop .Deop .st .PIE .dual .eppos .eppos2 .epneg

if nargin<3 || isempty(dual), dual = false; end
Top = PIE.T;    Aop = PIE.A;

% The inequality route declares no Deop slack, so there is nothing for the
% operator gate to verify against; refuse rather than return an uncheckable H.
% (build_stab_2d_st2 refuses the 2-D equivalent for the same reason.)
if isfield(st,'sosineq_on') && st.sosineq_on
    error('build_stab_1d:sosineq', ...
          'sosineq_on=1 has no Deop slack; the operator gate needs the equality route.');
end

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);

% ---- positive Lyapunov operator Pop (PDEstability STEP 1) ----------------
[prog, P1op] = poslpivar(prog, Top.dim, st.dd1, st.options1);
if st.override1~=1
    [prog, P2op] = poslpivar(prog, Top.dim, st.dd12, st.options12);
    Pop = P1op + P2op;
else
    Pop = P1op;
end
% strict positivity margin; eppos on the real-valued states, eppos2 on the
% distributed ones.  This is what floors max|Pop| and so keeps the gate's
% denominator from collapsing -- see pielr_opcheck.
Imat = blkdiag(st.eppos*eye(Pop.dim(1,:)),st.eppos2*eye(Pop.dim(2,:)));
Pop = Pop + mat2opvar(Imat, Pop.dim(:,2), PIE.vars, PIE.dom);

% ---- negativity operator Dop (PDEstability / _dual STEP 2) ---------------
if ~dual
    Dop = Top'*Pop*Aop + Aop'*Pop*Top + st.epneg*Top'*Pop*Top;
else
    Dop = Top*Pop*Aop' + Aop*Pop*Top' + st.epneg*Top*Pop*Top';
end

% ---- negativity slack Deop and the equality (PDEstability STEP 3) --------
[prog, De1op] = poslpivar(prog, Dop.dim, st.dd2, st.options2);
if st.override2~=1
    [prog, De2op] = poslpivar(prog, Dop.dim, st.dd3, st.options3);
    Deop = De1op + De2op;
else
    Deop = De1op;
end
prog = lpi_eq(prog, Dop+Deop, 'symmetric');

H.Top = Top;    H.Aop = Aop;    H.Pop = Pop;
H.Dop = Dop;    H.Deop = Deop;
H.eppos = st.eppos;   H.eppos2 = st.eppos2;   H.epneg = st.epneg;
H.dual = dual;  H.st = st;      H.PIE = PIE;
end

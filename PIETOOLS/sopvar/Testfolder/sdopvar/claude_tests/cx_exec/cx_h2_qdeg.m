function Qdeg = cx_h2_qdeg(PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QDEG = CX_H2_QDEG(PIE,ST) is the 1x3 degree array the 1-D non-coercive H2
% executives hand to lpivar for Q (PIETOOLS_H2_norm_c.m:123-134,
% PIETOOLS_H2_norm_o.m:120-131): get_lpivar_degs(Rop,Top) of the stock Rop
% declared at dd1/options1 (+ dd12/options12 when override1 ~= 1).
%
% Why: get_lpivar_degs accepts opvar/dopvar only (get_lpivar_degs.m:60,
% 114-115), so the container has no counterpart. Rop is declared on a
% scratch program and discarded; only its degrees are read, as
% test_lpivar_cdopvar_stability.m:108,133 does. lpivar_cdopvar's legacy
% [d1 d2 d3] mapping then reproduces lpivar's family.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Rop] = poslpivar(prog,PIE.T.dim,st.dd1,st.options1);
if st.override1~=1
    [~,P2op] = poslpivar(prog,PIE.T.dim,st.dd12,st.options12);
    Rop = Rop + P2op;
end
Qdeg = get_lpivar_degs(Rop,PIE.T);
end

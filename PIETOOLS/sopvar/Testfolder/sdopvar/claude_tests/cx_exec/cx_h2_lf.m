function [prog,P] = cx_h2_lf(prog,dm,sp,dom,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = CX_H2_LF(PROG,DM,SP,DOM,ST) declares the positive gramian /
% storage operator of the 1-D H2 executives over the spaces SP:
%
%   [prog,P1op] = poslpivar(prog,dim,dd1,options1);
%   if override1~=1, [prog,P2op] = poslpivar(prog,dim,dd12,options12);
%                    P = P1op+P2op;  else P = P1op;  end
%
% (e.g. PIETOOLS_H2_norm_c.m:123-130, PIETOOLS_H2_control.m:142-148).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[prog,P] = cx_h2_pos(prog,dm,sp,dom,st.dd1,st.options1);
if st.override1~=1
    [prog,P2] = cx_h2_pos(prog,dm,sp,dom,st.dd12,st.options12);
    P = P + P2;
end
end

function [prog,N] = cx_h2_slack(prog,D,dom,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,N] = CX_H2_SLACK(PROG,D,DOM,ST) declares the positive slack of the
% 1-D executives' equality branch (sosineq_on = 0) over the spaces of the
% square container D:
%
%   [prog,De1op] = poslpivar(prog,D.dim,dd2,options2);
%   if override2~=1, [prog,De2op] = poslpivar(prog,D.dim,dd3,options3);
%                    Deop = De1op+De2op;  else Deop = De1op;  end
%
% (e.g. PIETOOLS_H2_norm_c.m:166-177, PIETOOLS_H2_control.m:178-189). The
% spaces are read off D, not merged: poscopvar's single Gram over separate
% spaces gives the cone poslpivar gives over the merged R and L2 blocks
% (equal decision counts, test_copvar_kyp).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sp,dm] = cx_space_list(D,'out');
[prog,N] = cx_h2_pos(prog,dm,sp,dom,st.dd2,st.options2);
if st.override2~=1
    [prog,N2] = cx_h2_pos(prog,dm,sp,dom,st.dd3,st.options3);
    N = N + N2;
end
end

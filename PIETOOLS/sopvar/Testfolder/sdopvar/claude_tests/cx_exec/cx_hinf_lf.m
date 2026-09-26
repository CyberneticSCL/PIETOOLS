function [prog,Pm] = cx_hinf_lf(prog,Tm,PIE,st)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,PM] = CX_HINF_LF(PROG,TM,PIE,ST) declares the Lyapunov / storage
% operator of the 1-D H-infinity executives with poscopvar over the state
% spaces of TM (opvar2copvar(PIE.T)):
%
%   [prog,P1op] = poslpivar(prog,Top.dim,dd1,options1);
%   if override1~=1
%       [prog,P2op] = poslpivar(prog,Top.dim,dd12,options12);  P = P1+P2;
%   end
%
% (PIETOOLS_Hinf_gain.m:145-152, _coercive:133-140, _dual:149-156,
% _dual_coercive:148-155, Hinf_control:147-154, Hinf_estimator:146-153).
% Hinf_control passes Top.dim(:,1) for the psatz term, the same square
% spaces. No eppos here: the gain executives never add it, and the
% synthesis ones add it themselves after this call.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sp,dm] = cx_space_list(Tm,'out');      % zero-dimension spaces already dropped
[deg,co] = cx_hinf_posdeg(st.dd1,st.options1,sp);
[prog,Pm] = poscopvar(prog,dm,sp,PIE.dom,deg,co);
if st.override1~=1
    [deg,co] = cx_hinf_posdeg(st.dd12,st.options12,sp);
    [prog,P2] = poscopvar(prog,dm,sp,PIE.dom,deg,co);
    Pm = Pm + P2;
end
end

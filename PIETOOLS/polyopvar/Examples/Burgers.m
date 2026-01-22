%
% Burgers eq.
%
% PDE: u_t = u_ss + r*u - u*u_s;
% BCs: u(t,0) = u(t,1) = 0;
%
% fundamental state: v = u_ss \in L_2[0,1];
%
% maps from fundamental state to PDE states
% u   = (Tv) = \int_0^s T_1(s,t)v(t)dt + \int_s^1 T_2(s,t)v(t)dt;
% u_s = (Rv) = \int_0^s R_1(s,t)v(t)dt + \int_s^1 R_2(s,t)v(t)dt;
% T_1(s,t)=(s-1)*t; T_2(s,t) = s*(t-1);
% R_1(s,t)=t;       R_2(s,t) = t-1;
%
% PIE: (Tv_t) = v + (r*Tv) - (Tv)(Rv);
%
% PIE as a polyopvar: Z_1(v) = v; C = C_1 = [[C_111, C_112] 
%                                            [C_121, C_122]
%                                            [C_131, C_132]]
%
% C_111 = 1; C_112 = 1; C_121 = 1; C_122 = r*T; C_131 = -T; C_132 = R
% where nz = 1; m_1 = 3; d_1 = d_11 = 2; size(C_1) = (3,2)

% declare spatial variables, domain, and r.
pvar s t
dom = [0,1];
r = 1.0;

% Declare T as an opvar.
opvar T;
T1 = (s-1)*t; T.R.R1 = T1;
T2 = s*(t-1); T.R.R2 = T2;
T.I = dom;
T.var1 = s; T.var2 = t;

% Declare R as an opvar.
opvar R;
R1 = t;   R.R.R1 = R1;
R2 = t-1; R.R.R2 = R2;
R.I = dom;
R.var1 = s; R.var2 = t;

% Convert T, R to nopvar operators.
T = dopvar2ndopvar(T);
R = dopvar2ndopvar(R);

% Declare identity nopvar for identity operators present in PIE.
id = nopvar();
id.C = {1, 0, 0};
id.deg = 0;
id.dom = dom;
id.vars = [s, t];

% Declare PIE as polyopvar.
C1 = {id,  id;
      id, r*T;
      -T,   R};
C = tensopvar();
C.ops = {C1};
p = polyopvar('v',s,dom);
p.C = C;
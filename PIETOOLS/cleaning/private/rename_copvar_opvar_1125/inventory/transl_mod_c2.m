function T = transl_mod_c2(T,I)
% Probe copy of transl_mod: the check admits today's container name (as isa(T,'opvar') will after
% the rename), and each kernel is translated in a local K and assigned to the field once, so the
% container setter never receives a kernel in the auxiliary variables s_bar, theta_bar.
if nargin==1
    I = [-1,1];
elseif nargin>2
    error("Incorrect number of inputs; only 2 inputs are allowed");
end
if ~isa(T,'copvar') && ~isstruct(T)
    error("First argument must be a struct of kernels or an opvar object");
end
if length(I)~=2
    error("Second argument must be an array of length 2");
end

I_init = T.I;
a = I_init(1);
b = I_init(2);

c = I(1);
d = I(2);

s = T.var1;     theta = T.var2;
sbar = pvar([s.varname{1},'_bar']);
thetabar = pvar([theta.varname{1},'_bar']);

K = ((b-a)/(d-c))*subs(T.Q1, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
T.Q1 = subs(K, sbar, s);

K = subs(T.Q2, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
T.Q2 = subs(K, sbar, s);

K = subs(T.R.R0, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c));
T.R.R0 = subs(K, sbar, s);

K = ((b-a)/(d-c))*subs(subs(T.R.R1, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c)), theta, ((b-a)/(d-c))*thetabar+(a*d-b*c)/(d-c));
T.R.R1 = subs(subs(K, sbar, s),thetabar, theta);

K = ((b-a)/(d-c))*subs(subs(T.R.R2, s, ((b-a)/(d-c))*sbar+(a*d-b*c)/(d-c)), theta, ((b-a)/(d-c))*thetabar+(a*d-b*c)/(d-c));
T.R.R2 = subs(subs(K, sbar, s),thetabar,theta);

end

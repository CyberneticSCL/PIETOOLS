clc; clear;
pvar s t theta;
x = state('pde'); z = state('out'); w = state('in');

% heat equation
eqns = [diff(x,t)==diff(x,s,2)+(s^2-s)*w; 
             z==int(x,s,[0,1]);
             subs(x,s,0)==0; subs(x,s,1)==0];
% % transport equation
% eqns = [diff(x,t)==diff(x,s)+w; 
%              z==int(x,s,[0,1]);
%              subs(x,s,1)==0;];
% % reaction-diffusion equation
% eqns = [diff(x,t)==2*x+diff(x,s,2)+s*w; 
%              z==int(x,s,[0,1]);
%              subs(x,s,0)==0; subs(diff(x,s),s,1)==0];
%%
pde = sys();
pde = addequation(pde,eqns);
pie = convert(pde);
PIE = pie.params;
%%
st = lpisettings('heavy');
st.sos_opts.solver = 'mosek';
[prog,P,gam_hinf] = lpisolve(PIE,st,'l2gain');
[prog,P,gam_hinf_d] = lpisolve(PIE,st,'l2gain-dual');
[prog,P,gam_o] = lpisolve(PIE,st,'h2-o');
[prog,P,gam_c] = lpisolve(PIE,st,'h2-c');
%% find h2-norm & hinf-norm of the transfer function
lims = 200; dw = 0.05;
N = 2*lims/dw;
dW = linspace(-lims,lims,N);
intS = 0;
hinf_norm = 0;
for i=dW
    [tfncR,tfncI] = tf(PIE,i);
    H = double(tfncR.P)+j*double(tfncI.P);
    hinf_norm = max(hinf_norm,abs(double(tfncR.P)));
    tr = trace(H'*H);
    intS = intS+tr*dw;
end
h2norm = sqrt(intS/(2*pi));



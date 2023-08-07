clc; clear;
pvar s t theta;
x = state('pde'); z = state('out'); w = state('in');

% heat equation
eqns = [diff(x,t)==diff(x,s,2)+(s^2-s)*w; 
             z==int(x,s,[0,1]);
             subs(x,s,0)==0; subs(diff(x,s),s,1)==0];
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
opvar I; I.dim = PIE.T.dim; I.var1 = PIE.T.var1; I.var2 = PIE.T.var2;
I.R.R0 = 1;
lims = 50; dw = 0.5;
N = lims/dw;
dW = linspace(0.1,lims,N);
intS = 0;
hinf_norm = 0;
for i=1:length(dW)
    w = dW(i);
%     [tfncR,tfncI] = tf(PIE,w);
    % finding -A+i wT, analytically
    M_r = -inv(PIE.A+w^2*PIE.T*PIE.A*PIE.T);
    M_i = w*PIE.T*M_r;
    tfnewR = PIE.C1*M_r*PIE.B1;
    tfnewI = PIE.C1*M_i*PIE.B1;
    H = double(tfnewR.P)+1i*double(tfnewI.P);
    hinf_norm = max(hinf_norm,abs(H));
%     errR(i)=double(tfnewR.P-tfncR.P);
%     errI(i)=double(tfnewI.P-tfncI.P);
end
%%
plot(dW,errR,'x'); hold on;
plot(dW,errI,'o');
title('Error in approximation of transfer function.')
legend('Real part','Imaginary part');
xlabel('Frequency, w');
ylabel('err');


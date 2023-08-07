clc; clear;

pvar s t theta;
x = state('pde'); z = state('out'); w = state('in');

% heat equation
eqns = [diff(x,t)==diff(x,s,2); 
             z==int(x,s,[0,1]);
             subs(x,s,0)==0; subs(diff(x,s),s,1)==w];
% % transport equation
% eqns = [diff(x,t)==diff(x,s); 
%              z==int(x,s,[0,1]);
%              subs(x,s,1)==w;];
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
%% find h2-norm & hinf-norm of the transfer function
opvar I; I.dim = PIE.T.dim; I.var1 = PIE.T.var1; I.var2 = PIE.T.var2;
I.R.R0 = 1;
lims = 10; N=200;
dw = linspace(-lims,lims,N);
dW = 10.^(dw);
hinf_norm = 0;
for i=1:length(dW)
    w = dW(i);
    [tfncR,tfncI] = tf_2(PIE,w);
    H = double(tfncR)+1i*double(tfncI);
%     H = double(tfncR.P)+1i*double(tfncI.P);
    resp_num(i) = abs(H);
end
hinf_norm = max(resp_num);


loglog(dW,double(abs(resp_num)));
title('Bode plot: $\dot{\mathbf{x}}=\mathbf{x}_{ss},~~\mathbf{x}(0)=0,~~\mathbf{x}(1)=w,~~z=\int_0^1 \mathbf{x}(s) ds$','Interpreter','latex');
% title('Bode plot: $\dot{\mathbf{x}}=\mathbf{x}_{s},~~\mathbf{x}(1)=w,~~z=\int_0^1 \mathbf{x}(s) ds$','Interpreter','latex');
xlabel('$\omega$','Interpreter','latex');
ylabel('$| G(i \omega) |$','Interpreter','latex');
str = "Max |G| = "+num2str(hinf_norm);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
%%
syms sx;
tf_analytic = (cosh(sqrt(sx))-1)/(sqrt(sx)*sinh(sqrt(sx)));
figure(2);
w = linspace(-10,10,500);
freq = 1i*10.^(w);
resp = subs(tf_analytic,sx,freq);
hinf_analytic = max(double(abs(resp)));
loglog(-1i*freq,double(abs(resp)));
title('Bode plot: $\dot{\mathbf{x}}=\mathbf{x}_{ss},~~\mathbf{x}(0)=0,~~\mathbf{x}(1)=w,~~z=\int_0^1 \mathbf{x}(s) ds$','Interpreter','latex');
xlabel('$\omega$','Interpreter','latex');
ylabel('$| G(i \omega) |$','Interpreter','latex');
str = "Max |G| = "+num2str(hinf_analytic);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
%%
syms sx;
tf_analytic = (1-exp(-sx))/sx;
figure(3);
w = linspace(-10,10,500);
freq = 1i*10.^(w);
resp = subs(tf_analytic,sx,freq);
hinf_analytic = max(double(abs(resp)));
loglog(-1i*freq,double(abs(resp)));
title('Bode plot: $\dot{\mathbf{x}}=\mathbf{x}_{s},~~\mathbf{x}(1)=w,~~z=\int_0^1 \mathbf{x}(s) ds$','Interpreter','latex');
xlabel('$\omega$','Interpreter','latex');
ylabel('$| G(i \omega) |$','Interpreter','latex');
str = "Max |G| = "+num2str(hinf_analytic);
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
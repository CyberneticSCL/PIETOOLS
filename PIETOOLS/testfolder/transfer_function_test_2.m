clc; clear;

pvar s t theta;
x = state('pde'); z = state('out'); w = state('in');

% heat equation
eqns = [diff(x,t)==diff(x,s,2)+w; 
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
st = lpisettings('heavy');
st.sos_opts.solver = 'mosek';
[prog,P,gam_hinf] = lpisolve(PIE,st,'l2gain');
%%
clear resp_num;
lims = 3; Ns=50; N = 10;
dw = linspace(-lims,lims,N);
dW = 10.^(dw);
% dW = [-dW(end:-1:1),dW];
hinf_norm = 0;
opvar Z Zb;
Z.dim = PIE.C1.dim;
Zb.dim = PIE.B1.dim;
C = PIE.C1; B = PIE.B1; Tw = PIE.Tw; Ap = PIE.A; T = PIE.T;
nz = PIE.C1.dim(1,1);
nw = PIE.B1.dim(2,1);
si = linspace(C.I(1),C.I(2),Ns);
for i=1:length(dW)
    i
    w = dW(i);
    invop = inv_opvar_new([-Ap,w*T; -w*T,-Ap],1e-8,Ns);
    Cq = [C,Z;Z,C]; Bq = [B-w*Tw,Zb;Zb,B-w*Tw];
    [tfncR,tfncI] = findfinalTF(Cq.Q1,invop,Bq.Q2,si);
%     bigTF = Cq*invop*Bq;
%     tmp = bigTF.P;
%     tfncR = tmp(1:nz,1:nw); tfncI =tmp(1:nz,nw+1:2*nw);
%     [tfncR,tfncI] = tf_2(PIE,w);
    H = double(tfncR)+1i*double(tfncI);
    resp_num(i) = abs(H);
end
hinf_norm = max(resp_num);

% syms sx;
% tf_analytic = (cosh(sqrt(sx))-1)/(sqrt(sx)*sinh(sqrt(sx)));
% resp = subs(tf_analytic,sx,1i*dW);
% hinf_analytic = max(double(abs(resp)));
% loglog(dW,double(abs(resp))); hold on;
plot(dW,double(abs(resp_num)));
title('$\dot{\mathbf{x}}=\mathbf{x}_{s}+w,~~\mathbf{x}(1)=0,~~z=\int_0^1 \mathbf{x}(s) ds$','Interpreter','latex');
xlabel('$\omega$','Interpreter','latex');
ylabel('$| G(i \omega) |$','Interpreter','latex');
% legend('Analytical','Numerical');
%%
tmp =invop{2};
tmp = reshape(cat(3,tmp{:}),[2,2,Ns,Ns]);

subplot(2,2,1);
surf(squeeze(tmp(1,1,:,:))); hold on;
subplot(2,2,2);
surf(squeeze(tmp(1,2,:,:))); 
subplot(2,2,3);
surf(squeeze(tmp(2,1,:,:))); 
subplot(2,2,4);
surf(squeeze(tmp(2,2,:,:))); 
colorbar;

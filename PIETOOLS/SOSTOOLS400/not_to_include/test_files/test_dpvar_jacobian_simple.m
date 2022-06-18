clc; clear; 
% Scalable Global Stability Analysis Test Code


% pvar p1 p2
% g=(1-p1^2-p2^2); % norm{p1,p2} \le 1
% % g1=p1*(1-p2);
% % g2=p1*(1-p1)*p2*(1-p2); % p1,p2 \in [0,1]

n=5;deg=4;
xvartable=[];
for i=1:n
c{i}=['p',int2str(i)];
temp=pvar(c{i});
xvartable=[xvartable temp];
end
Z=monomials(xvartable,1:deg);
rng(5);
mons=randi(length(Z),n,1);
Zs=Z(mons);
f_new=-jacobian(Zs,xvartable)*xvartable';

%f_temp=-xvartable';

disp('setting up the SOS program using pvar + finding jacobian/lyapunov derivative')
tic
prog1 = sosprogram(xvartable);
[prog1,V]=sossosvar(prog1,Z);
Vd=jacobian(V,xvartable)*f_new;
toc



disp('setting up the SOS program using dpvar + finding jacobian/lyapunov derivative')
tic
prog2 = sosprogram(xvartable);
[prog2,V2]=DPsossosvar(prog2,Z);
Vd2=DPjacobian(V2,xvartable)*f_new;
toc

disp('Testing equality of lyapunov derivative using dpvar2poly and poly minus');
tic;
tmp = dpvar2poly(Vd2)-Vd;
toc;
max(max(abs(tmp.coefficient)))

% disp('Testing equality of lyapunov derivative using poly2dpvar and dpvar minus');
% tic;
% tmp2 = Vd2 - poly2dpvar(Vd);
% toc;

max(max(abs(tmp2.C)))


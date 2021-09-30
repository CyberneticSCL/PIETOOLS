clc; clear; 
% setup"
% We have a linear system of the form
% \dot x(t)=A(p)x(t) 
% where we would like to prove stability for all p \in G:=\{p: g(p)\ge 0\}
% Construct P(p) > \epsilon I           for all p \in G
% such that A(p)^T P(p)+P(p)A(p)<0      for all p \in G
%

pvar p1 p2
g=(1-p1^2-p2^2); % norm{p1,p2} \le 1
% g1=p1*(1-p2);
% g2=p1*(1-p1)*p2*(1-p2); % p1,p2 \in [0,1]

n=10;
deg=2;
ep=.01;

A=-eye(n)+.25*tril(ones(n))*p1-.25*triu(ones(n))*p2;

mastervartable=[p1 p2];

Z=monomials_nd(mastervartable,0:deg,n);
Z2=monomials_nd(mastervartable,0:(deg+1),n);
nZ1=size(Z,1);
nZ2=size(Z2,1);

prog1 = sosprogram(mastervartable);
disp('Robust Stability Analysis using dpvar using sosposmatr')
disp('creating Positive matrix variables')
tic
[prog1,P1a] = DPsosposmatr(prog1,nZ1);
[prog1,P1b] = DPsosposmatr(prog1,nZ1);
[prog1,P1c] = DPsosposmatr(prog1,nZ2);

Pp1a=Z'*P1a*Z+ep*eye(n);
Pp1b=Z'*P1b*Z;
Pp1c=Z2'*P1c*Z2;

disp('Constructing derivative matrix')
Ineq=A'*Pp1a+Pp1a*A;
% Ineq=-Pp1b-Pp1c*g
temp=Ineq+Pp1b+Pp1c*g;
[prog1]=DPsoseq(prog1, temp);% DPsosmateq does not exist

disp('Solving SOS program')
%prog1_s=sossolve(prog1)
toc 


% disp('Robust Stability Analysis using pvar using sosposmatr')
% disp('creating Positive matrix variables')
% tic
% prog2 = sosprogram(mastervartable);
% [prog2,P2a] = sosposmatr(prog2,nZ1);
% [prog2,P2b] = sosposmatr(prog2,nZ1);
% [prog2,P2c] = sosposmatr(prog2,nZ2);
% 
% Pp2a=Z'*P2a*Z+ep*eye(n); % currently bugged
% Pp2b=Z'*P2b*Z;
% Pp2c=Z2'*P2c*Z2;
% 
% disp('Constructing derivative matrix')
% Ineq=A'*Pp2a+Pp2a*A;
% % Ineq=-Pp2b-Pp2c*g
% temp=Ineq+Pp2b+Pp2c*g;
% % tic 
%  [prog2_tmp]=sosmateq(prog2, temp); 
% % toc
% %tic
% %[prog2]=soseq(prog2, temp);  %I'm not sure this actually works for matrices...

disp('Solving SOS program')

%prog2_s=sossolve(prog2)
toc

disp('Robust Stability Analysis using dpvar using sospolymatrixvar')

tic
Zsym=monomials(mastervartable,0:2*deg);
prog3 = DPsosprogram(mastervartable);
[prog3,Pp3a] = DPsospolymatrixvar(prog3,Zsym,[n n],'symmetric');
[prog3] = DPsosmatrixineq(prog3,Pp3a-ep*eye(n));
[prog3,Pp3c] = DPsospolymatrixvar(prog3,Zsym,[n n],'symmetric');
[prog3] = DPsosmatrixineq(prog3,Pp3c);
Ineq=A'*Pp3a+Pp3a*A;
%Ineq+Pp3c*g < 0
[prog3] = DPsosmatrixineq(prog3,-Pp3c*g-Ineq);
toc
prog3temp=DPsossolve(prog3)

disp('Robust Stability Analysis using dpvar using sosposmatrvar')
% strangely, this version has almost twice as many variables.
%Zsym=monomials(mastervartable,0:2*deg);
tic
prog3 = DPsosprogram(mastervartable);
[prog3,Pp3a] = DPsosposmatrvar(prog3,n,2*deg,mastervartable);
Pp3a=Pp3a+ep*eye(n);
%[prog3] = DPsosmatrixineq(prog3,Pp3a-ep*eye(n));
[prog3,Pp3c] = DPsosposmatrvar(prog3,n,2*deg,mastervartable);
%[prog3] = DPsosmatrixineq(prog3,Pp3c);
Ineq=A'*Pp3a+Pp3a*A;
%Ineq+Pp3c*g < 0
[prog3] = DPsosmatrixineq(prog3,-Pp3c*g-Ineq);
toc
prog3temp=DPsossolve(prog3)

disp('Robust Stability Analysis using pvar using sospolymatrixvar')

tic
Zsym=monomials(mastervartable,0:2*deg);
prog3 = sosprogram(mastervartable);
[prog3,Pp3a] = sospolymatrixvar(prog3,Zsym,[n n],'symmetric');
[prog3] = sosmatrixineq(prog3,Pp3a-ep*eye(n));
[prog3,Pp3c] = sospolymatrixvar(prog3,Zsym,[n n],'symmetric');
[prog3] = sosmatrixineq(prog3,Pp3c);
Ineq=A'*Pp3a+Pp3a*A;
%Ineq+Pp3c*g < 0
[prog3] = sosmatrixineq(prog3,-Pp3c*g-Ineq);
toc
%prog3temp=sossolve(prog3)

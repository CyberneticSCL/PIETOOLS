clc; clear; 

pvar th ksi
pvar gam
n=10;
deg=4;
g=th*(1-th);
%g=th*(1-th)*ksi*(1-ksi);

mastervartable=[th ksi];

disp('Error testing and timing for declaration of Positive Matrix Variable')
Z=monomials_nd(mastervartable,0:deg,n);
nZ1=size(Z,1);

prog1 = sosprogram(mastervartable,[gam]);
prog2 = sosprogram(mastervartable,[gam]);
prog3 = sosprogram(mastervartable);
disp('creating Positive matrix as dpvar')
tic
[prog1,P1] = DPsosposmatr(prog1,nZ1);
[prog1,P1_s] = DPsosposmatr_struct(prog1,6, eye(6));
toc
disp('creating Positive matrix as pvar')
tic
[prog2,P2] = sosposmatr(prog2,nZ1);
toc

disp('Constructing SOS variable using sossosmatrixvar - may take a while')

Zsym=monomials(mastervartable,2*deg);
tic
[prog3,P] = sospolymatrixvar(prog3,Zsym,[n n],'symmetric');
[prog3] = sosmatrixineq(prog3,P);
toc
disp('Error testing and timing for construction of Positive Matrix with variables')
disp('testing that these positive matrices are the same')
tic
% P1_dpvar=poly2dpvar(P1);
% toc
% tic
% error=P1-P2;
% toc
% error.C


disp('Constructing SOS variable as dpvar')

tic
%[prog1,P1] = DPsosposmatr(prog1,nZ1);
%Mvar1L=Z'*P1;
%Mvar1R=P1*Z;
Mvar1=Z'*P1*Z;
%Mvar1=Mvar1L*Z;
toc
disp('Constructing SOS variable as pvar')

tic
%Mvar2L=Z'*P2;
%Mvar2R=P2*Z;
%[prog2,P2] = sosposmatr(prog2,nZ1);
Mvar2=Z'*P2*Z;
%Mvar2=Mvar2L*Z;
toc
%Mvar2L
%error4=Mvar2L-Mvar1L;
%error4=Mvar2L-Mvar1L;

% Mvar1_pvar=dpvar2poly(Mvar1);
% error_pvar=Mvar1_pvar-Mvar2;
% error_pvar.C

% error_dpvar=Mvar1-Mvar2;
% error_dpvar.C

% error2=P1-P2;
% error2.C


disp('Multiplying SOS variable by g as dpvar')

tic
Mvar1g=Mvar1*g;
%Mvar1g_alt=Mvar1*(eye(n)*g);
toc

% Mvar1g_alt=mtimes_edit(Mvar1,g);

disp('Multiplying SOS variable by g as pvar')

tic
Mvar2g=Mvar2*g;
toc

%err_alt=Mvar1g_alt-Mvar2g;

% Mvar1g_pvar=dpvar2poly(Mvar1g);
% 
% errorMg_pvar=Mvar1g_pvar-Mvar2g;
% 
% errorMg=Mvar2g-Mvar1g;






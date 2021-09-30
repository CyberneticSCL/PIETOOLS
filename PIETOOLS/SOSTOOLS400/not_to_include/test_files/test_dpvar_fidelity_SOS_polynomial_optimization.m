clc; clear; 
% Scalable Global Polynomial Optimization Fidelity Test Code


% % % Unconstrained optimization problem %  % 
% f(x) is the Goldstein-Price function

pvar x1 x2
f1 = x1 + x2 +1;
f2 = 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2;
f3 = 2*x1 - 3*x2;
f4 = 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2;
f = (1 + f1^2*f2)*(30 + f3^2*f4);

% Solution: fmin = 3
%   x = [0, -1];

vartable = [x1 x2];

disp(' - Solving the unconstrained optimization problem using pvar - ')
disp('setting up the sos and poly variables')
tic
prog1 = sosprogram(vartable);
tp1 = toc;
disp(['   elapsed time is ',num2str(tp1),'s'])

disp('setting up the decision variable gamma')
pvar gam;                       % let gamma be a lower bound on the value of f
tic
prog1 = sosdecvar(prog1, gam);  % set gam = gamma as decision var
prog1 = sossetobj(prog1, -gam); % minimize -gamma
tp2 = toc;
disp(['   elapsed time is ',num2str(tp2),'s'])

disp('setting up the inequality constraint f>=gamma')
expr = (f-gam);
tic
prog1 = sosineq(prog1,expr);
tp3 = toc;
disp(['   elapsed time is ',num2str(tp3),'s'])

disp('solving the SOS program using pvar')
tic
prog1_s = sossolve(prog1);
tp4 = toc;
disp(['   elapsed time is ',num2str(tp4),'s'])
fminp0 = sosgetsol(prog1_s,gam);

% % %

disp(' ')
disp(' - Solving the unconstrained optimization problem using dpvar - ')
disp('setting up the sos and poly variables')
tic
prog2 = sosprogram(vartable);
td1 = toc;
disp(['   elapsed time is ',num2str(td1),'s'])

disp('setting up the decision variable gamma')
pvar gam;                       % let gamma be a lower bound on the value of f
tic
prog2 = DPsosdecvar(prog2, gam);  % set gam = gamma as decision var
prog2 = DPsossetobj(prog2, -gam); % minimize -gamma
td2 = toc;
disp(['   elapsed time is ',num2str(td2),'s'])

disp('setting up the inequality constraint f>=gamma')
expr = (f-gam);
tic
prog2 = DPsosineq(prog2,expr);
td3 = toc;
disp(['   elapsed time is ',num2str(td3),'s'])

disp('solving the SOS program using dpvar')
tic
prog2_s = DPsossolve(prog2);
td4 = toc;
disp(['   elapsed time is ',num2str(td4),'s'])
fmind0 = DPsosgetsol(prog2_s,gam);


tptot = tp1 + tp2 + tp3;% + tp4;
tdtot = td1 + td2 + td3;% + td4;

disp(' ')
if double(fminp0) == double(fmind0)
    disp(['The poly and DP implementations produce the same minimal value gamma = ',num2str(double(fmind0))])
end
%disp(['Set-up for unconstrained problem using DP works ',num2str(tptot/tdtot),' times faster than poly version']);
%disp(['Solving the unconstrained problem using DP works ',num2str(tp4/td4),' times faster than poly version']);


% % % Simple constrained optimization problem % % %

% min   f(x)
% s.t.  g(x)>=0
pvar x y

f = x^4 + y^4 - 2*y*x^3 - 3*y^2*x^2 + 150*(x^2 + y^2);
g1 = 12^2 - x^2;
g2 = 12^2 - y^2;
g3 = 2*12^2 - (x^2 + y^2);

% Solution: fmin = -19008
%   x = [-12, -12];

deg = 3;
vartable = [x y];
Z = monomials(vartable,0:deg);


disp(' ')
disp(' - Solving the simple optimization problem using pvar - ')
disp('setting up the sos and poly variables')
tic
prog1 = sosprogram(vartable);
[prog1,s1] = sossosvar(prog1,Z);
[prog1,s2] = sossosvar(prog1,Z);
[prog1,s3] = sossosvar(prog1,Z);
tp1 = toc;
disp(['   elapsed time is ',num2str(tp1),'s'])

disp('setting up the decision variable gamma')
pvar gam;                       % let gamma be a lower bound on the value of f
tic
prog1 = sosdecvar(prog1, gam);  % set gam = gamma as decision var
prog1 = sossetobj(prog1, -gam); % minimize -gamma
tp2 = toc;
disp(['   elapsed time is ',num2str(tp2),'s'])

disp('setting up the inequality constraint f>=gamma')
expr = (f-gam) - s1*g1 - s2*g2 - s3*g3;
tic
prog1 = sosineq(prog1,expr);
tp3 = toc;
disp(['   elapsed time is ',num2str(tp3),'s'])

disp('solving the SOS program using pvar')
tic
prog1_s = sossolve(prog1);
tp4 = toc;
disp(['   elapsed time is ',num2str(tp4),'s'])
fminp1 = sosgetsol(prog1_s,gam);

% % %

disp(' ')
disp(' - Solving the simple optimization problem using dpvar - ')
disp('setting up the sos and poly variables')
tic
prog2 = DPsosprogram(vartable);
[prog2,s1] = DPsossosvar(prog2,Z);
[prog2,s2] = DPsossosvar(prog2,Z);
[prog2,s3] = DPsossosvar(prog2,Z);
td1 = toc;
disp(['elapsed time is ',num2str(td1),'s'])

disp('setting up the decision variable gamma')
dpvar gam;                       % let gamma be a lower bound on the value of f
tic
prog2 = DPsosdecvar(prog2, gam);  % set gam = gamma as decision var
prog2 = DPsossetobj(prog2, -gam); % minimize -gamma
td2 = toc;
disp(['   elapsed time is ',num2str(td2),'s'])

disp('setting up the inequality constraint f>=gamma')
expr = (f-gam) - s1*g1 - s2*g2 - s3*g3;
tic
prog2 = DPsosineq(prog2,expr);
td3 = toc;
disp(['   elapsed time is ',num2str(td3),'s'])

disp('solving the SOS program using dpvar')
tic
prog2_s = DPsossolve(prog2);
td4 = toc;
disp(['   elapsed time is ',num2str(td4),'s'])
fmind1 = DPsosgetsol(prog2_s,gam);


tptot = tp1 + tp2 + tp3;% + tp4;
tdtot = td1 + td2 + td3;% + td4;

disp(' ')
if double(fminp1) == double(fmind1)
    disp(['The poly and DP implementations produce the same lower bound gamma = ',num2str(double(fmind1))])
end
%disp(['Set-up for simple constrained problem using DP works ',num2str(tptot/tdtot),' times faster than poly version']);
%disp(['Solving the simple constrained problem using DP works ',num2str(tp4/td4),' times faster than poly version']);



% % % Difficult problem % % %

% Example from Hesameddin Mohammadi and Matthew Peet, 2017:
% min   f(x)
% s.t.  g(x)>=0,    h(x)=0

pvar x1 x2 x3 x4 x5 x6

f = 7*x1*x5^3 + 6*x1*x5^2*x6 + 9*x2*x4^3 + 4*x2*x4*x5 + 3*x2*x5*x6 + x3*x4*x5;

g1 = 100 - (x1^2 + x2^2 + x3^2 + x4^2 + x5^2 + x6^2);
g2 = x1^3 + x2^2*x4 + x3*x5^2;
g3 = x2^2*x1 + x5^3 + x4*x1*x2;
h1 = x1 + x2^2 - x3^2 + x4*x5;
h2 = x5*x1 - x4^2;

% Solution: fmin = -3693.3
%   x = [5.1416, 3.9307, 0.7568, −4.6777, 4.2676, −4.1504]


%deg_f = 2;
deg_g1 = 2;     deg_g2 = 2;     deg_g3 = 2;
deg_h1 = 4;     deg_h2 = 4;
vartable = [x1 x2 x3 x4 x5 x6];
%Z_f = monomials(vartable,0:deg_f);
Z_g1 = monomials(vartable,0:deg_g1);
Z_g2 = monomials(vartable,0:deg_g2);
Z_g3 = monomials(vartable,0:deg_g3);
Z_h1 = monomials(vartable,0:deg_h1);
Z_h2 = monomials(vartable,0:deg_h2);

disp('  ')
disp(' - Solving the more difficult optimization problem using pvar - ')
disp('setting up the sos and poly variables')
tic
prog1 = sosprogram(vartable);
[prog1,s1] = sossosvar(prog1,Z_g1);
[prog1,s2] = sossosvar(prog1,Z_g2);
[prog1,s3] = sossosvar(prog1,Z_g3);
[prog1,t1] = sospolyvar(prog1,Z_h1);
[prog1,t2] = sospolyvar(prog1,Z_h2);
tp1 = toc;
disp(['   elapsed time is ',num2str(tp1),'s']) 

disp('setting up the decision variable gamma')
pvar gam;                       % let gamma be a lower bound on the value of f
tic
prog1 = sosdecvar(prog1, gam);  % set gam = gamma as decision var
prog1 = sossetobj(prog1, -gam); % minimize -gamma
tp2 = toc;
disp(['   elapsed time is ',num2str(tp2),'s']) 
% gam = -3700;
% tp2 = 0;

disp('setting up the inequality constraint f>=gamma')
tic
expr = (f-gam) - s1*g1 - s2*g2 - s3*g3 - t1*h1 - t2*h2;
toc
tic
prog1 = sosineq(prog1,expr);
tp3 = toc;
disp(['   elapsed time is ',num2str(tp3),'s']) 

disp('solving the SOS program using pvar')
tic
prog1_s = sossolve(prog1);
tp4 = toc;
disp(['   elapsed time is ',num2str(tp4),'s']) 
fminp2 = sosgetsol(prog1_s,gam);

% % %

disp(' ')
disp(' - Solving the more difficult optimization problem using dpvar - ')
disp('setting up the sos and poly variables')
tic
prog2 = DPsosprogram(vartable);
[prog2,s1] = DPsossosvar(prog2,Z_g1);
[prog2,s2] = DPsossosvar(prog2,Z_g2);
[prog2,s3] = DPsossosvar(prog2,Z_g3);
[prog2,t1] = DPsospolyvar(prog2,Z_h1);
[prog2,t2] = DPsospolyvar(prog2,Z_h2);
td1 = toc;
disp(['   elapsed time is ',num2str(td1),'s'])

disp('setting up the decision variable gamma')
dpvar gam;                       % let gamma be a lower bound on the value of f
tic
prog2 = DPsosdecvar(prog2, gam);  % set gam = gamma as decision var
prog2 = DPsossetobj(prog2, -gam); % minimize -gamma
td2 = toc;
disp(['   elapsed time is ',num2str(td2),'s'])
% gam = -3700;
% td2 = 0;

disp('setting up the inequality constraint f>=gamma')
tic
expr = (f-gam) - s1*g1 - s2*g2 - s3*g3 - t1*h1 - t2*h2;
toc
tic
prog2 = DPsosineq(prog2,expr);
td3 = toc;
disp(['   elapsed time is ',num2str(td3),'s'])

disp('solving the SOS program using dpvar')
tic
prog2_s = DPsossolve(prog2);
td4 = toc;
disp(['   elapsed time is ',num2str(td4),'s'])
fmind2 = DPsosgetsol(prog2_s,gam);



tptot = tp1 + tp2 + tp3;% + tp4;
tdtot = td1 + td2 + td3;% + td4;

disp(' ')
if double(fminp2) == double(fmind2)
    disp(['The poly and DP implementations produce the same lower bound gamma = ',num2str(double(fmind2))])
end
%disp(['Set-up for more difficult problem using DP works ',num2str(tptot/tdtot),' times faster than poly version']);
%disp(['Solving the more difficult problem using DP works ',num2str(tp4/td4),' times faster than poly version']);
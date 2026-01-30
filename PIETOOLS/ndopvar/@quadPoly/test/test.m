% Example: build a 2x2 quadPoly with small monomial bases and print it
clc; clear;

nvars = {'s1','s2'}; 
nvart = {'t1','t2'};
dim = [2,2];     % matrix dimension
nmons = [5,5];  % number of monomials in s, theta
maxdeg = [3,3]; % maximum degree in s, theta
den = 0.05;     % sparsity in coefficients
F = quadPoly.randquadPoly(dim,nmons,nvars,nvart,maxdeg,den);
G = quadPoly.randquadPoly(dim,nmons,nvars,nvart,maxdeg,den);


%% addition
clc;
H1 = F+G;
H2 = G-G;

disp(H1);
disp(H2);
%% transpose
clc;
disp(F);
disp(F');

%% mtimes
clc;
disp(F);
disp(F');
disp(F'*F);


%% concatenation
clc;
disp([F, F]);
disp([F;F]);


%% indexing
clc;
disp(F);
disp(F(4));
disp(F(2:end,1:end));
disp(F(end));

%% 
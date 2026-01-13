clc; clear;
n = 4;

dim = [1,1];
vars_in=[];
for i=1:n
    eval("pvar s"+num2str(i));
    vars_in = [vars_in, eval("s"+num2str(i))];
end
vars_out = vars_in(1);
%%

mons = monomials(vars_in,0:5);

k = (0:3^n-1)';          % 3^n-by-1
pow3 = 3.^(0:n-1);       % 1-by-n
keysallidx = mod(floor(k ./ pow3), 3);   % 3^n-by-n, digits in base-3

keysall = nDopvar.index2keys(keysallidx,n);

params = cell(1,length(keysall));
for i=1:length(keysall)
    coef = rand(size(mons,1),prod(dim));
    params{i} = polynomial(coef,mons.degmat,mons.varname, dim);
end
paramsDict = dictionary(keysall,params);
P = nDopvar(vars_in.varname, vars_out.varname, dim, paramsDict);

%%
% addition 
disp('Timing addition');
fadd = @() plus(P,P);
tadd = timeit(fadd);
fprintf('Average addition time for %d D to %d D PI operator: %d\n', n,n,tadd);


% transpose
disp('Timing Adjoint');
fadj = @() ctranspose(P);
tadj = timeit(fadj);
fprintf('Average Adjoint time for %d D to %d D PI operator: %d\n', n,n,tadj);

R= P';

% composition
disp('Timing Composition 1');
fmul = @() mtimes(R,P);
tmul = timeit(fmul);
fprintf('Average Composition time (method 1) for %d D to %d D PI operator: %d\n', n,n,tmul);

disp('Timing Composition 2');
fmul2 = @() mtimes2(R,P);
tmul2 = timeit(fmul2);
fprintf('Average Composition time (method 2) for %d D to %d D PI operator: %d\n', n,n,tmul2);


% Test_cg_sopvar_nd - tests cg_sopvar on sopvar operators in n = 1,2,3 spatial

%
% Cost grows steeply with n (about 40 s at n=3 for 8 iterations).
%
% Initial coding AT - 10/13/2026

maxit = 10;
tol_ip = 1e-8;      % relative tolerance on the inner products
tol_err = 1e-8;     % tolerance on ||X-Xt||^2
m = 2;

rng(1); 

for N = [2]
    % Declare N variables and their domain
    % N = 2 is 5-10 mins
    % N = 3 -- out of memory
    d = 2;
    var1_name = [repmat('s',[N,1]),num2str((1:N)')];
    var2_name = [var1_name,repmat('_dum',[N,1])];
    var1 = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
    var2 = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));
    dom = [zeros(N,1),ones(N,1)];
    
    A_nopvar = rand_ndopvar([m,m],d,dom,var1,var2,cell(0, 1));
    A = nopvar2sopvar(A_nopvar);
    % X_nopvar = rand_ndopvar([m,0],d,dom,var1,var2,cell(0, 1));
    % X_nopvar = rand_ndopvar([m,0],d,dom,var1,var2,cell(0, 1));

    % random sopvar mx1 with 0 input variables
    vars = struct();  vars.in  = {};     vars.out = cellstr(var1_name).';     % no input variables
    domB  = struct();  domB.in   = zeros(0,2);  domB.out = dom;
    ZL   = A.ZL;          % monomial degrees in s1 and s2 (column vectors)
    nZ   = prod(cellfun(@numel,ZL));  % 3*2 = 6 monomials
    C    = {sparse(randn(m*nZ,1))};     % one parameter, nZ x 1, since there is no input
    X   = sopvar(C,vars,ZL,{},domB,[m,1]);   % ZR is empty
 
    B = A*X; 

    I1= id_sopvar(m,A.vars,A.dom);
    A = A'*A + I1; 
    B = A'*B;
    [Xr,res,it] = cg_sopvar_pd(A,B,1e-10,40); % min ||A X - B ||

    % [Xr,res,it] = cg_sopvar(A,B,1e-10,40);
    E = A*Xr-B;   W = E'*E;   err = full(W.params{1});
    fprintf('n=%d passed: it=%d relRes=%.1e err=%.1e\n',N,it,res(end)/res(1),err);
end

 
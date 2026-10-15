% kernel_degree_test.m -- test kernel_degree(F) and kernel_degree2(F).
%
% For the degree d and spatial degree V_mon of SOS_DP, the spatial degree
% formula used in kernel_degree is evaluated directly from the degrees of
% the kernels of DP, and compared with kernel_degree and kernel_degree2.

% DJ, 10/15/2026: Updated to use new polyopvar functions;

clear stateNameGenerator;
pvar s t
dom = [0,1];
u   = pde_var(s,dom);
PDE = [diff(u,t)==diff(u,s,2) + 5*u + u^2;
       subs(u,s,dom(1))==0;
       subs(u,s,dom(2))==0];
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
x   = PIE.x;

% each row: [d,opdeg], the degree of the basis and the spatial degree in SOS_DP
test_cases = [1,0; 1,1; 1,2; 2,1; 2,2];
nfail = 0;
for k=1:size(test_cases,1)
    d = test_cases(k,1);    V_mon = test_cases(k,2);
    prog = piesos_program(x);
    [prog,DP] = SOS_DP(prog,d,V_mon,x,dom);

    % The formula used in kernel_degree, with the degrees of the kernels
    nz = size(DP.degmat,1); % number of operators.
    deg_opts = zeros(1,nz); % max degrees of each operator.
    for i=1:nz
        Kc = kernels(DP.C.ops{1,i});
        for l=1:numel(Kc)
            if (isa(Kc{l},'polynomial') || isa(Kc{l},'dpvar')) && ~isempty(Kc{l}.degmat)
                % max. of sum over rows of degmat = degree of polynomial.
                deg_opts(i) = max(deg_opts(i),max(sum(full(Kc{l}.degmat),2)));
            end
        end
    end
    [max_deg, arg_deg] = max(deg_opts);
    dd = full(DP.degmat(arg_deg)); % degree of FDP in linear form - should always be even for SOS and LFs.
    dd = ceil(dd/2); % degree of FDP in quadratic form - used in formula.
    opdeg = (max_deg-1) / (2*dd); % this could be negative if deg_opts = 0
    opdeg = max(0, opdeg); % maximum polynomial degree.
    opdeg = ceil(opdeg);

    od1 = kernel_degree(DP);
    od2 = kernel_degree2(DP);
    % The formula and the functions agree; SOS_DP uses monomials of degree
    % opdeg in s and in s_dum separately (kron), and kernel_degree returns 2*opdeg, not opdeg
    % (the kernels of DP have about twice the degree assumed in the formula, see kernel_degree2)
    ok =(od1==2*V_mon) && (od2==od1) && (opdeg==od1);
    fprintf('d=%d, V_mon=%d: formula %d, kernel_degree %d, kernel_degree2 %d: %s\n',d,V_mon,opdeg,od1,od2,string(ok));
    nfail = nfail + ~ok;
end
if nfail==0
    fprintf(' --- kernel_degree tests passed ---\n');
else
    error('kernel_degree_test:Failed','%d of the cases did not return V_mon.',nfail);
end

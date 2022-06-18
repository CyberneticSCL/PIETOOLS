% Code testing proper functioning of DPposlpivar
% The code compares the results to those attain building the positive 
% dopvar directly using quadratic products
% A comparison with the poly poslpivar would maybe make more sense, but
% some situations seem not to be supported by this function

% NOTE: Test does not work now that we are using fastint2str

tol = 1e-15;
ntests = 50;

% maximal degree of the monomials
d0_max = 2;     d1_max = 3;     d2_max = 4;
% maximal size of the dopvar
n1_max = 5;     n2_max = 5;
% options in calling the function
test_p = 3;     test_s = 3;     test_e = 5;
% NOTE: diags option is not supported at this point...

if test_p==0
    psatz = 0;
elseif test_p==1
    psatz = 1;
end

if test_s==0
    sep=0;
elseif test_s==1
    sep=1;
end

if test_e==0
    excludeL = [0,0,0,0];
elseif test_e==1
    excludeL = [1,0,0,0];
elseif test_e==2
    excludeL = [0,1,0,0];
elseif test_e==3
    excludeL = [0,0,1,0];
elseif test_e==4
    excludeL = [0,0,0,1];
end

% % % Run a number of (random) tests % % %

for testnumber = 1:ntests

if test_p>=3
    psatz = randi(2)-1;
end
if test_s>=3
    sep = randi(2)-1;
end
if test_e>=5
    excludeL = randi(2,1,4) - 1;
    if all(excludeL(1:4-sep))
        excludeL(randi(4-sep)) = 0;
    end
end    

% % % Set up a random example % % %

I = sort(2*rand(1,2) - 1,'ascend');    % domain of the variables s,theta

n1 = randi(n1_max);
n2 = randi(n2_max);
n = [n1,n2];        % size of problem

d0 = randi(d0_max);
d1 = randi(d1_max,1,3);     d1(end) = max(d1);
d2 = randi(d2_max,1,3);     d2(end) = max(d2);
d = {d0,d1,d2};     % degrees of monomials Z1s, Z2sth, Z3sth

options.psatz = psatz;          % toggle to use psatz or not
options.exclude = excludeL;     % binary array indicating monomials to exclude
options.sep = sep;              % toggle to set R2=R1

pvar s theta
var1 = s;   var2 = theta;
prog = sosprogram([s,theta]);


% % % Build the poslpivars % % %

[~,P1,LLL] = DPposlpivar_directv1(prog,n,I,d,options);
[~,P2] = DPposlpivar(prog,n,I,d,options);   %<- no subsref version
%[~,P3] = DPposlpivar_nosubsref(prog,n,I,d,options);
[~,P4] = DPposlpivarQ(prog,n,I,d,options);  %<- quadvar version



% % % Check correctness of the results % % %

if all(excludeL(1:3)==[0,0,0]) && (all(d1==d2) || excludeL(4) == 1 || sep) % in these cases, we can compare with poly version
    [~,P0] = poslpivar(prog,n,I,d,options);

elseif 0 % ~psatz
% If psatz=0, we can directly construct the positive dopvar using the
% monomials and LLL
% doesn't work yet....

    Zvar = opvar();

    if ~excludeL(1)
        Zvar.P = eye(n(1));
    end

    bZ1s=[];
    if ~excludeL(2)
        % Constructing Z1(s)
        nZ1=d{1}+1;
        Z1degmat = (0:d{1})';
        Z1coeff = speye(nZ1);
        Z1varname = var1.varname;
        Z1matdim = [nZ1 1];
        Z1s=polynomial(Z1coeff,Z1degmat,Z1varname,Z1matdim);

        nBZ1=n2*nZ1;

        for i=1:n2
            bZ1s=blkdiag(bZ1s,Z1s);
        end

        Zvar.R.R0 = bZ1s;
        Zvar.R.R1 = zeros(size(bZ1s));
        Zvar.R.R2 = zeros(size(bZ1s));
    else
        nZ1 = 0;
    end

    bZ2sth=[];
    if ~excludeL(3)
    % Constructing Z2(s,th)
    % In this implementation, Z2 will have degree d{2,2} in th and degree d{2,1} in s
    % and max degree of s+th is d{2,3}. Similarly for Z3(s,th)

        Z2degmat = [repmat((0:d{2}(1))',d{2}(2)+1,1),vec(repmat(0:d{2}(2),d{2}(1)+1,1))];
        Z2degmat(sum(Z2degmat,2)>d{2}(3),:)= [];
        nZ2 = size(Z2degmat,1);
        Z2coeff = speye(nZ2);
        Z2varname = [var1.varname; var2.varname];
        Z2matdim = [nZ2 1];
        Z2sth=polynomial(Z2coeff,Z2degmat,Z2varname,Z2matdim);

        for i=1:n2
            bZ2sth=blkdiag(bZ2sth,Z2sth);
        end

        Zvar.R.R0 = [Zvar.R.R0; zeros(size(bZ2sth))];
        Zvar.R.R1 = [Zvar.R.R1; bZ2sth];
        Zvar.R.R2 = [Zvar.R.R2; zeros(size(bZ2sth))];
        
    else
        nZ2 = 0;
    end

    bZ3sth=[];
    nZ3 = 0;
    if sep
        Zvar.R.R0 = [Zvar.R.R0; zeros(size(bZ2sth))];
        Zvar.R.R1 = [Zvar.R.R1; zeros(size(bZ2sth))];
        Zvar.R.R2 = [Zvar.R.R2; bZ2sth];
    elseif ~excludeL(4)
    % Constructing Z3(s,th)    
        Z3degmat = [repmat((0:d{3}(1))',d{3}(2)+1,1),vec(repmat(0:d{3}(2),d{3}(1)+1,1))];
        Z3degmat(sum(Z3degmat,2)>d{3}(3),:)= [];
        nZ3=size(Z3degmat,1);
        Z3coeff = speye(nZ3);
        Z3varname = [var1.varname; var2.varname];
        Z3matdim = [nZ3 1];
        Z3sth=polynomial(Z3coeff,Z3degmat,Z3varname,Z3matdim);

        for i=1:n2
            bZ3sth=blkdiag(bZ3sth,Z3sth);
        end

        Zvar.R.R0 = [Zvar.R.R0; zeros(size(bZ3sth))];
        Zvar.R.R1 = [Zvar.R.R1; zeros(size(bZ3sth))];
        Zvar.R.R2 = [Zvar.R.R2; bZ3sth];
    end

    nBZ1=n2*nZ1;
    nBZ2=n2*nZ2;
    nBZ3=n2*nZ3;
    dimL1=(1-excludeL(1))*n1;
    dimL2=(1-excludeL(2))*nBZ1;
    dimL3=(1-excludeL(3))*nBZ2;
    dimL4=(1-excludeL(4))*nBZ3;
    dimLLL=dimL1+dimL2+dimL3+dimL4;
    ind{1}=1:dimL1; 
    ind{2}=(dimL1+1):(dimL1+dimL2); 
    ind{3}=(dimL1+dimL2+1):(dimL1+dimL2+dimL3);
    ind{4}=(dimL1+dimL2+dimL3+1):dimLLL;

    if sep
        LLL(ind{1},ind{4}) = LLL(ind{1},ind{3});
        LLL(ind{2},ind{4}) = LLL(ind{2},ind{3});
        LLL(ind{3},ind{4}) = LLL(ind{3},ind{3});
        LLL(ind{4},ind{4}) = LLL(ind{3},ind{3});
        LLL(ind{4},ind{1}) = LLL(ind{3},ind{1});
        LLL(ind{4},ind{2}) = LLL(ind{3},ind{2});
        LLL(ind{4},ind{3}) = LLL(ind{3},ind{3});
    end


    P0 = Zvar' * LLL * Zvar;

else
    % with no better alternative, use the direct version for comparison
    P0 = P1;
    
end
    
% Check correctness of standard DPposlpivar result
difP_2 = P0.P - P2.P;
difQ1_2 = P0.Q1 - P2.Q1;
difQ2_2 = P0.Q2 - P2.Q2;
difR0_2 = P0.R.R0 - P2.R.R0;
difR1_2 = P0.R.R1 - P2.R.R1;
difR2_2 = P0.R.R2 - P2.R.R2;

if max(difP_2.C,[],'all') >= tol
    error('DPposlpivar does not produce the correct P matrix')
end
if max(difQ1_2.C,[],'all') >= tol || max(difQ2_2.C,[],'all') >= tol
    error('DPposlpivar does not produce the correct Q1 or Q2')
end
if max(difR0_2.C,[],'all') >= tol
    error('DPposlpivar does not produce the correct R0')
end
if max(difR1_2.C,[],'all') >= tol || max(difR2_2.C,[],'all') >= tol
    error('DPposlpivar does not produce the correct R1 or R2')
end


% % Check correctness of DPposlpivar without subsref result
% difP_3 = P0.P - P3.P;
% difQ1_3 = P0.Q1 - P3.Q1;
% difQ2_3 = P0.Q2 - P3.Q2;
% difR0_3 = P0.R.R0 - P3.R.R0;
% difR1_3 = P0.R.R1 - P3.R.R1;
% difR2_3 = P0.R.R2 - P3.R.R2;
% 
% if max(difP_3.C,[],'all') >= tol
%     error('DPposlpivar_nosubsref does not produce the correct P matrix')
% end
% if max(difQ1_3.C,[],'all') >= tol || max(difQ2_3.C,[],'all') >= tol
%     error('DPposlpivar_nosubsref does not produce the correct Q1 or Q2')
% end
% if max(difR0_3.C,[],'all') >= tol
%     error('DPposlpivar_nosubsref does not produce the correct R0')
% end
% if max(difR1_3.C,[],'all') >= tol || max(difR2_3.C,[],'all') >= tol
%     error('DPposlpivar_nosubsref does not produce the correct R1 or R2')
% end

% Check correctness of DPposlpivar with quadvar result

if 0 %any(excludeL)
    P_4 = P4.P;         P_1 = P1.P;
    Q1_4 = P4.Q1;       Q1_1 = P1.Q1;
    Q2_4 = P4.Q2;       Q2_1 = P1.Q2;
    R0_4 = P4.R.R0;     R0_1 = P1.R.R0;
    R1_4 = P4.R.R1;     R1_1 = P1.R.R1;
    R2_4 = P4.R.R2;     R2_1 = P1.R.R2;
    
    P_4.dvarname = P_1.dvarname;
    Q1_4.dvarname = Q1_1.dvarname;
    Q2_4.dvarname = Q2_1.dvarname;
    R0_4.dvarname = R0_1.dvarname;
    R1_4.dvarname = R1_1.dvarname;
    R2_4.dvarname = R2_1.dvarname;
    
    difP_4 = P_1 - P_4;
    difQ1_4 = Q1_1 - Q1_4;
    difQ2_4 = Q2_1 - Q2_4;
    difR0_4 = R0_1 - R0_4;
    difR1_4 = R1_1 - R1_4;
    difR2_4 = R2_1 - R2_4;
    
else   
    
    difP_4 = P0.P - P4.P;
    difQ1_4 = P0.Q1 - P4.Q1;
    difQ2_4 = P0.Q2 - P4.Q2;
    difR0_4 = P0.R.R0 - P4.R.R0;
    difR1_4 = P0.R.R1 - P4.R.R1;
    difR2_4 = P0.R.R2 - P4.R.R2;
    
end

if max(difP_4.C,[],'all') >= tol
    error('DPposlpivar_quadvar does not produce the correct P matrix')
end
if max(difQ1_4.C,[],'all') >= tol || max(difQ2_4.C,[],'all') >= tol
    error('DPposlpivar_quadvar does not produce the correct Q1 or Q2')
end
if max(difR0_4.C,[],'all') >= tol
    error('DPposlpivar_quadvar does not produce the correct R0')
end
if max(difR1_4.C,[],'all') >= tol || max(difR2_4.C,[],'all') >= tol
    error('DPposlpivar_quadvar does not produce the correct R1 or R2')
end



end

disp('No errors were encountered testing DPposlpivar!')


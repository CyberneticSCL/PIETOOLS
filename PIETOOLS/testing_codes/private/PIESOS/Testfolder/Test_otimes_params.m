% Test_otimes_params.m -- compare otimes_params(A,B,otype) with the
% product computed directly from the parameters A and B, for random
% parameters of random sizes. The product type otype and the type of the
% parameters clss ('double', 'sparse', 'poly', 'quadPoly', 'dpvar',
% 'dpvar2' or 'quadPolyCell') are set in the script; 'quadPolyCell' (cells
% of 'polynomial' parameters expanded in monomials by quad_expand) is the
% default. The result C and the products A and B are applied to random
% vectors, and the difference of the results is displayed as nrm; no
% pass/fail check is made.
%
% DJ (CC), 10/15/2026: Initial coding;

% Declare the type of product to compute
otype = [1,2];
m1 = randi(5);      n1 = randi(5);
if otype(1)==1
    m2 = m1;
else
    m2 = randi(5);
end
if otype(2)==1
    n2 = n1;
else
    n2 = randi(5);
end

% Set parameters
ndvars = 10;   % number of decision variables
M1 = randi([0,5]);        N1 = randi([0,5]);
M2 = randi([0,5]);        N2 = randi([0,5]);
deg1 = 4;
deg2 = 4;
nZ1 = 20;
nZ2 = 20;

% Set the variables
svars = cell(1,M1+M2);
for i=1:M1+M2
    svars{i} = ['s',num2str(i)];
end
tvars = cell(1,N1+N2);
for i=1:N1+N2
    tvars{i} = ['t',num2str(i)];
end
s_idcs1 = randperm(M1+M2,M1);       t_idcs1 = randperm(N1+N2,N1);
s_idcs2 = randperm(M1+M2,M2);       t_idcs2 = randperm(N1+N2,N2);
svar1 = svars(s_idcs1);             tvar1 = tvars(t_idcs1);
svar2 = svars(s_idcs2);             tvar2 = tvars(t_idcs2);

pdegs1 = randi([0,deg1],[N1+M1,1]);
pdegs2 = randi([0,deg2],[N2+M2,1]);

% Generate random arrays of desired dimensions
clss = 'quadPolyCell';
if strcmp(clss,'double')
    A = rand([m1,n1]);
    B = rand([m2,n2]);
elseif strcmp(clss,'sparse')
    A = sprand(m1,n1,0.4);
    B = sprand(m2,n2,0.6);
elseif strcmp(clss,'poly') || strcmp(clss,'quadPoly')
    A_poly = rand_poly([m1,n1],[svar1,tvar1]',pdegs1,nZ1);
    B_poly = rand_poly([m2,n2],[svar2,tvar2]',pdegs2,nZ2);
    if strcmp(clss,'quadPoly')
        A = quadPoly.polynomial2quadPoly(A_poly,svar1,tvar1);
        B = quadPoly.polynomial2quadPoly(B_poly,svar2,tvar2);
    else
        A = A_poly;
        B = B_poly;
    end
elseif strcmp(clss,'dpvar')
    A = rand_dpvar([m1,n1],ndvars,[svar1,tvar1],pdegs1,nZ1);
    B = rand_poly([m2,n2],[svar2,tvar2],pdegs2,nZ2);
elseif strcmp(clss,'dpvar2')
    A = rand_poly([m1,n1],N1,deg1,nZ1);
    B = rand_dpvar([m2,n2],ndvars,N2,deg2,nZ2);
elseif strcmp(clss,'quadPolyCell')
    A = cell(3,2);      A_poly = cell(3,2);
    B = cell(2,3);      B_poly = cell(2,3);
    for i=1:numel(A_poly)
        A_poly{i} = rand_poly([m1,n1],[svar1,tvar1]',pdegs1,nZ1);
        B_poly{i} = rand_poly([m2,n2],[svar2,tvar2]',pdegs2,nZ2);

        %A{i} = quadPoly.polynomial2quadPoly(A_poly{i},svar1,tvar1);
        %B{i} = quadPoly.polynomial2quadPoly(B_poly{i},svar2,tvar2);
    end
    [Acell,ZA] = quad_expand(A_poly,struct('out',{svar1},'in',{tvar1}));
    [Bcell,ZB] = quad_expand(B_poly,struct('out',{svar2},'in',{tvar2}));

    A = struct();           B = struct();
    A.C = Acell;            B.C = Bcell;
    A.Zs = ZA.out;          B.Zs = ZB.out;
    A.Zt = ZA.in;           B.Zt = ZB.in;
    A.ns = svar1;           B.ns = svar2;
    A.nt = tvar1;           B.nt = tvar2;
    A.dim = [m1,n1];        B.dim = [m2,n2];
end

% Take the product
C = otimes_params(A,B,otype);

if strcmp(clss,'quadPoly')
    A = A_poly;
    B = B_poly;
    C = quadPoly.quadPoly2polynomial(C);
end
if strcmp(clss,'quadPolyCell')
    A = A_poly;
    B = B_poly;
    Cstruct = C;
    C = C.C;
    for k=1:numel(C)
        C{k} = quadPoly(C{k},Cstruct.Zs,Cstruct.Zt,Cstruct.dim,Cstruct.ns,Cstruct.nt);
        C{k} = quadPoly.quadPoly2polynomial(C{k});
    end
end


% Generate random vectors
if otype(2)==2
    if ~isa(A,'cell')
        vA = rand([n1,1]);
        vB = rand([n2,1]);
        vC = kron(vA,vB);
        Av = A*vA;
        Bv = B*vB;
        Cv = C*vC;
    else
        vA = repmat({rand([n1,1])},size(A));
        vB = repmat({rand([n2,1])},size(B));
        vC = repmat({kron(vA{1},vB{1})},size(C));
        Av = 0;
        for i=1:numel(vA)
            %vA{i} = rand([n1,1]);
            Av = Av + A{i}*vA{i};
        end
        Bv = 0;
        for j=1:numel(vB)
            %vB{j} = rand([n2,1]);
            Bv = Bv + B{j}*vB{j};
        end
        Cv = 0;
        Cv_true = 0;
        K = ndims(A);
        szB = size(B,1:K);
        for k=1:numel(vC)
            % Determine which parameters in A and B correspond to parameter
            % k in C
            idcsC = cell(1,K);
            [idcsC{:}] = ind2sub(size(C,1:K),k);
            idcsC = cell2mat(idcsC);
            idcsA = ceil(idcsC./szB);
            idcsB = idcsC - (idcsA-1).*szB;
            idcsA = num2cell(idcsA);
            idcsB = num2cell(idcsB);
            % Take the Kronecker product of the appropriate terms
            %Ak = vA{idcsA{:}};
            %Bk = vB{idcsB{:}};
            %vC{k} = kron(Ak,Bk);
            % 
            Cv = Cv + C{k}*vC{k};
        end
    end
elseif otype(1)==2
    vA = rand([1,m1]);
    vB = rand([1,m2]);
    vC = kron(vA,vB);
    Av = vA*A;
    Bv = vB*B;
    Cv = vC*C;
end

% Compare the products
if all(otype==2)
    err = Cv - repelem(Av,m2,1).*repmat(Bv,m1,1);
elseif any(otype==2)
    err = Cv - Av.*Bv;
else
    err = C - A.*B;
end
if strcmp(clss,'double') || strcmp(clss,'sparse')
    nrm = norm(full(err))
else
    nrm = max(max(abs(err.C)))
end


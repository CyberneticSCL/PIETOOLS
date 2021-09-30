
clear
clc

% size of problem
n1 = 1;     n2 = 8;

% degrees of monomials
%d = {3,[3,4,5],[3,4,5]}; 
d = {2,[3,4,5],[3,4,5]}; 

% Domain
I = [0,1];


% Set up a program
pvar s theta sss
prog = sosprogram([s,theta]);

var1 = s;
var2 = theta;

% Construct monomial vector Z1(s)
nZ1=d{1}+1;
Z1degmat = (0:d{1})';
Z1coeff = speye(nZ1);
Z1varname = var1.varname;
Z1matdim = [nZ1 1];
Z1s=polynomial(Z1coeff,Z1degmat,Z1varname,Z1matdim);

% Construct monomial vector Z2(s,theta)
Z2degmat = [repmat((0:d{2}(1))',d{2}(2)+1,1),vec(repmat(0:d{2}(2),d{2}(1)+1,1))];
Z2degmat(sum(Z2degmat,2)>d{2}(3),:)= [];
nZ2=size(Z2degmat,1);
Z2coeff = speye(nZ2);
Z2varname = [var1.varname; var2.varname];
Z2matdim = [nZ2 1];
Z2sth=polynomial(Z2coeff,Z2degmat,Z2varname,Z2matdim);

nBZ1=n2*nZ1;
nBZ2=n2*nZ2;

% Create blockdiagonal versions of monomial vectors
bZ1s=[];
for i=1:n2
    bZ1s=blkdiag(bZ1s,Z1s);
end
bZ2sth=[];
for i=1:n2
    bZ2sth=blkdiag(bZ2sth,Z2sth);
end


bZ1th=subs(bZ1s,var1,var2);

bZ2ths=var_swap(bZ2sth,var1,var2);
bZ2etath=subs(bZ2sth,var1,sss);
bZ2seta = subs(bZ2sth,var2,sss);
bZ2etas=subs(bZ2ths,var2,sss);

% % % Set up dpvar implementation % % %

% Create the dpvar decision variable object LLL_dp
dimL1 = n1;
dimL2 = nBZ1;
dimL3 = nBZ2;
dimL4 = nBZ2;
dimLLL=dimL1+dimL2+dimL3+dimL4;
disp('Setting up the decision matrix as dpvar...')
tic
[~,LLL_dp]=DPsosposmatr(prog,dimLLL);
toc

% Split the object into the different subcomponents Q{i,j}
ind{1}=1:dimL1; 
ind{2}=(dimL1+1):(dimL1+dimL2); 
ind{3}=(dimL1+dimL2+1):(dimL1+dimL2+dimL3);
ind{4}=(dimL1+dimL2+dimL3+1):dimLLL;

disp(' ')
disp('Decomposing decision dpvar into the necessary blocks without compression...')
Q_dp = cell(4,4);
tic
for i=1:4
    for j=1:4
        if i>=j         
            Q_dp{i,j}=LLL_dp(ind{i},ind{j});
            if i~=j
                Q_dp{j,i}=Q_dp{i,j}.';
            end
        end
    end
end
toc

disp('Decomposing decision dpvar into the necessary blocks with compression...')
Q_dp_com = cell(4,4);
tic
for i=1:4
    for j=1:4
        if i>=j         
            Q_dp_com{i,j}=DPcompress(LLL_dp(ind{i},ind{j}));
            if i~=j
                Q_dp_com{j,i}=Q_dp_com{i,j}.';
            end
        end
    end
end
toc

% % % Set up polynomial implementation % % %

% Create the polynomial decision variable object LLL_poly
disp(' ')
disp('Setting up the decision matrix as polynomial...')
tic
[~,LLL_poly]=sosposmatr(prog,dimLLL);
toc

% Split the object into the different subcomponents Q{i,j}
ind{1}=1:dimL1; 
ind{2}=(dimL1+1):(dimL1+dimL2); 
ind{3}=(dimL1+dimL2+1):(dimL1+dimL2+dimL3);
ind{4}=(dimL1+dimL2+dimL3+1):dimLLL;

disp('Decomposing decision polynomial into the necessary blocks...')
Q_poly = cell(4,4);
tic
for i=1:4
    for j=1:4
        if i>=j         
            Q_poly{i,j}=LLL_poly(ind{i},ind{j});
            if i~=j
                Q_poly{j,i}=Q_poly{i,j}.';
            end
        end
    end
end
toc

%--------------------------------------------------------------------------

% Performing the necessary multiplication using dpvars and the polynomial
% subroutines
disp(' ')
disp('Computing products using dpvars...')
tic
Z1QZ2 = bZ1s'*Q_dp{2,3}*bZ2sth;
toc
tic
Z2QZ1 = bZ2ths'*Q_dp{4,2}*bZ1th;
toc
tic
Z2QZ2 = bZ2etas'*Q_dp{3,3}*bZ2etath;
toc

disp('Computing products using dpvar subroutines...')
tic
Z1QZ2_alt = Z1QZ2construct(Z1s,Z2sth,Q_dp{2,3},n2);
toc
tic
Z2QZ1_alt = Z1QZ2construct(Z2sth,Z1s,Q_dp{4,2},n2);
toc
tic
Z2QZ2_alt = Z2QZ3construct(Z2sth,Z2sth,Q_dp{3,3},n2,var1,sss);
toc

disp(' ')
disp('Computing products using compressed dpvars...')
tic
Z1QZ2_com = bZ1s'*Q_dp_com{2,3}*bZ2sth;
toc
tic
Z2QZ1_com = bZ2ths'*Q_dp_com{4,2}*bZ1th;
toc
tic
Z2QZ2_com = bZ2etas'*Q_dp_com{3,3}*bZ2etath;
toc

disp('Computing products using compressed dpvars with subroutines...')
tic
Z1QZ2_alt_com = Z1QZ2construct(Z1s,Z2sth,Q_dp_com{2,3},n2);
toc
tic
Z2QZ1_alt_com = Z1QZ2construct(Z2sth,Z1s,Q_dp_com{4,2},n2);
toc
tic
Z2QZ2_alt_com = Z2QZ3construct(Z2sth,Z2sth,Q_dp_com{3,3},n2,var1,sss);
toc

disp(' ')
disp('Computing the products using polynomial subroutines...')
tic
Z1QZ2_poly = constructA1(Q_poly{2,3},Z1s,Z2sth,1,n2,nZ1,nZ2,var1,var2);
toc
tic
Z2QZ1_poly = constructB1(Q_poly{4,2},Z1s,Z2sth,1,n2,nZ1,nZ2,var1,var2); 
toc
tic
Z2QZ2_poly = constructCA1(Q_poly{3,3},Z2sth,1,n2,nZ2,var1,var2,sss); 
toc

disp(' ')
% Integrating the dpvars and polynomials
disp('Performing the necessary integration using dpvars...')
tic
Z2QZ2_int = int(Z2QZ2,sss,var1,I(2));
toc
disp('Performing the necessary integration using compressed dpvars...')
tic
Z2QZ2_com_int = int(Z2QZ2_com,sss,var1,I(2));
toc
disp('Performing the necessary integration using polynomials...')
tic
Z2QZ2_poly_int = int(Z2QZ2_poly,sss,var1,I(2));
toc

disp(' ')
% Adding the different terms (having to establish common bases) takes a lot
% of time using dpvars...
disp('Adding the products using dpvars...')
tic
R1_dp = Z1QZ2;
R1_dp = R1_dp + Z2QZ1;
R1_dp = R1_dp + Z2QZ2_int;
toc

disp('Adding the products using compressed dpvars...')
tic
R1_dp_com = Z1QZ2;
R1_dp_com = R1_dp_com + Z2QZ1;
R1_dp_com = R1_dp_com + Z2QZ2_int;
toc

disp('Adding the products using polynomials...')
tic
R1_poly = Z1QZ2_poly;
R1_poly = R1_poly + Z2QZ1_poly; 
R1_poly = R1_poly + Z2QZ2_poly_int; 
toc


%error1 = dpvar2poly(Z1QZ2) - Z1QZ2_poly;
%error2 = dpvar2poly(Z2QZ1) - Z2QZ1_poly;
%error3 = dpvar2poly(Z2QZ2) - Z2QZ2_poly;

%tol = 1e-15;
%error0 = dpvar2poly(R1_dp) - R1_poly;
%if max(error0.coeff,[],'all')>tol
%    error('R1_dp and R1_poly do not match!')
%end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % Polynomial subroutines % % %

function [A1] = constructA1(Q12,ZZZth,ZZZthksi,g,n,nZth,nZthksi,var1,var2)

avarname=[Q12.varname; ZZZthksi.varname ];

if strcmp(var1.varname{1},ZZZthksi.varname{1}) % var1 is in the first position
        bZdegmat1=repmat([ZZZth.degmat zeros(nZth,1)],n,1);
else
        bZdegmat1=repmat([zeros(nZth,1) ZZZth.degmat],n,1);
end
    bZdegmat2=repmat(ZZZthksi.degmat,n,1);

[PIlist, PJlist]=find(Q12.coefficient); % This is the list of terms in P, but PIlist and Pjlist will be relatively short. 
[Prowu,Pcolu] = ind2sub([n*nZth n*nZthksi],PJlist); % this returns the matrix locations for every term in P
adegmat = [Q12.degmat(PIlist,:) bZdegmat1(Prowu,:)+bZdegmat2(Pcolu,:)];
  

[row,col] = ind2sub([n*nZth n*nZthksi],PJlist);
newrow=ceil(row/nZth);
newcol=ceil(col/nZthksi);
newidx=sub2ind([n n],newrow, newcol);

nPIlist=(1:length(PIlist))';
coeff=sparse(nPIlist,newidx,1);
amatdim=[n n];
A1=g*polynomial(coeff,adegmat,avarname,amatdim);
end

function [B1] = constructB1(Q31,ZZZth,ZZZthksi,g,n,nZth,nZthksi,var1,var2)
tvarname=[Q31.varname; ZZZthksi.varname{2}; ZZZthksi.varname{1}];
bZdegmat1=repmat(ZZZthksi.degmat,n,1);
bZdegmat2=repmat([ZZZth.degmat zeros(nZth,1)],n,1);

nZ1=nZthksi;nbZ1=n*nZ1;
nZ2=nZth;nbZ2=n*nZ2;

[PIlist, PJlist]=find(Q31.coefficient);  
[Prowu,Pcolu] = ind2sub([nbZ1 nbZ2],PJlist); 
tdegmat = [Q31.degmat(PIlist,:) bZdegmat1(Prowu,:)+bZdegmat2(Pcolu,:)];
[row,col] = ind2sub([nbZ1 nbZ2],PJlist);
newrow=ceil(row/nZ1);
newcol=ceil(col/nZ2);
newidx=sub2ind([n n],newrow, newcol);
nPIlist=(1:length(PIlist))';
coeff=sparse(nPIlist,newidx,1);
tmatdim=[n n];

B1=g*polynomial(coeff,tdegmat,tvarname,tmatdim);
end

function CA1 = constructCA1(Q33,ZZZthksi,g,n,nZthksi,var1,var2,dum)
tvarname=[Q33.varname; var1.varname{1}; var2.varname{1}; dum.varname{1}];


bZdegmat1=repmat([ZZZthksi.degmat(:,2) zeros(nZthksi,1)  ZZZthksi.degmat(:,1)],n,1); %bZomth
bZdegmat2=repmat([zeros(nZthksi,1) ZZZthksi.degmat(:,2)  ZZZthksi.degmat(:,1)],n,1); %bZomksi
nZ1=nZthksi;nbZ1=n*nZ1;
nZ2=nZthksi;nbZ2=n*nZ2;

[PIlist, PJlist]=find(Q33.coefficient); % This is the list of terms in P, but PIlist and Pjlist will be relatively short. 
[Prowu,Pcolu] = ind2sub([nbZ1 nbZ2],PJlist); % this returns the matrix locations for every term in P
tdegmat = [Q33.degmat(PIlist,:) bZdegmat1(Prowu,:)+bZdegmat2(Pcolu,:)];


newrow=ceil(Prowu/nZ1);
newcol=ceil(Pcolu/nZ2);

newidx=sub2ind([n n],newrow, newcol);
nPIlist=(1:length(PIlist))';
coeff=sparse(nPIlist,newidx,1);
tmatdim=[n n];
CA1=g*polynomial(coeff,tdegmat,tvarname,tmatdim);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % dpvar subroutines % % %


function [ZQZ] = Z1QZ1construct(Z1s,Q,n)
% This function computes the product bZ1s'*Q*bZ1s, producing a dpvar of
% dimensions [n,n];

nc = size(Q,2);
Zdeg = Z1s.degmat;
Zntt = size(Zdeg,1);

% We first build the new degmat
[row,~] = find(triu(repmat((1:Zntt)',1,Zntt))); % row indices of elements of Q
[~,col] = find(triu(repmat((1:Zntt),Zntt,1)));  % column indices
degmat = Zdeg(row,:) + Zdeg(col,:); % degmat of the product Z1s'*Q*Z1s
ntt = size(degmat,1);   nd = length(Q.dvarname)+1;

% Extract the coefficients and their locations from Q (Cval will be just
% ones)
[Ci,~,Cval] = find(Q.C);

% Each column in Q.C must be shifted to multiply with the appropriate row
% in the degmat
jndx_mat = 1 + cumsum(0:1:(Zntt-1)) + (0:1:(Zntt-1))';  % matrix of row indices in degmat coefficients Q(1:Zntt,1:Zntt) should multiply with
jndx_mat = triu(jndx_mat) + triu(jndx_mat,1)';          % account for the symmetry in Q
jndx_mat = repmat(jndx_mat,[n,n]);                      % matrix of indices indicating which row in degmat each coefficient Q(1:end,1:end) should be multiplied with
jndx_mat = jndx_mat + kron((0:ntt:(n-1)*ntt),ones(1,Zntt)); % adjust according to column position in the final object
Cjnew = vec(jndx_mat);                                  % adjusted column indices for coefficients in C

% Different coefficients Q(i,j) will add to the same element of the final
% object, we have to compress the row-indices accordingly
i_cmprss = repmat((0:nd:nd*(Zntt-1))',n,1) + kron((0:nd*(Zntt-1):(n-1)*nd*(Zntt-1))',ones(Zntt,1));
Cinew = Ci - repmat(i_cmprss,nc,1);

% Now we can build the new coefficient matrix, and corresponding dpvar
Cnew = sparse(Cinew,Cjnew,Cval,nd*n,ntt*n);
ZQZ = dpvar(Cnew,degmat,Z1s.varname,Q.dvarname,[n,n]);

end



function [ZQZ] = Z1QZ2construct(ZL,ZR,Q,n)
% This function computes the product bZLs'*Q*bZRsth or bZLths'*Q*bZRth, 
% producing a dpvar of dimensions [n,n];
% Input Z1 must be a monomial on just s or s and theta
% Input Z2 must be a monomial on s and theta or just s

% can probably still be more efficient

if length(ZL.varname)==1
    var1 = ZL.varname;      % 's'
    varname = ZR.varname;   % {'s';'theta'}
else
    var1 = ZR.varname;      % 's'
    varname = ZL.varname;   % {'s';'theta'}
end

nc = size(Q,2);
ZLdeg = ZL.degmat;
ZLntt = size(ZLdeg,1);
ZRdeg = ZR.degmat;
ZRntt = size(ZRdeg,1);

% We first build the new degmat
%row = repmat((1:ZLntt)',ZRntt,1);           % row indices of elements of Q
%col = kron((1:ZRntt)',ones(ZLntt,1));       % column indices

ntt = ZLntt * ZRntt;        nd = length(Q.dvarname)+1;
% degmat of the product Z1s'*Q*Z2sth
if length(ZL.varname)==1
    if strcmp(var1{1},ZR.varname{1}) % var1 (s) is in the first position
        %degmat = [ZLdeg(row,:),zeros(length(row),1)] + ZRdeg(col,:);
        degmat = [repmat(ZLdeg,ZRntt,1),zeros(ntt,1)] + ...
                    kron(ZRdeg,ones(ZLntt,1));
    else
        %degmat = [zeros(length(row),1),ZLdeg(row,:)] + ZRdeg(col,:);
        degmat = [zeros(ntt,1),repmat(ZLdeg,ZRntt,1)] + ...
                    kron(ZRdeg,ones(ZLntt,1));
    end
else
    if strcmp(var1{1},ZL.varname{1}) % var1 (s) is in the first position
        %degmat = ZLdeg(row,2:-1:1) + [zeros(length(col),1),ZRdeg(col,:)];
        degmat = repmat(ZLdeg(:,2:-1:1),ZRntt,1) + ...
                    [zeros(ntt,1),kron(ZRdeg,ones(ZLntt,1))];
    else
        %degmat = ZLdeg(row,:) + [ZRdeg(col,:),zeros(length(col),1)];
        degmat = repmat(ZLdeg,ZRntt,1) + ...
                    [kron(ZRdeg,ones(ZLntt,1)),zeros(ntt,1)];
    end
end
%ntt = size(degmat,1);      nd = length(Q.dvarname)+1;


% Extract the coefficients and their locations from Q (Cval will be just
% ones)
[Ci,~,Cval] = find(Q.C);

% Each column in Q.C must be shifted to multiply with the appropriate row
% in the degmat
jndx_mat = reshape((1:ntt*n)',ZLntt,n*ZRntt);   % matrix of row indices in repeated degmat coefficients Q(1:ZLntt,1:end) should multiply with
jndx_mat = repmat(jndx_mat,[n,1]);              % matrix of indices indicating which row in degmat each coefficient Q(1:end,1:end) should be multiplied with
Cjnew = vec(jndx_mat);                          % adjusted column indices for coefficients in C

% Different coefficients Q(i,j) will add to the same element of the final
% object, we have to compress the row-indices accordingly
i_cmprss = repmat((0:nd:nd*(ZLntt-1))',n,1) + kron((0:nd*(ZLntt-1):(n-1)*nd*(ZLntt-1))',ones(ZLntt,1));
Cinew = Ci - repmat(i_cmprss,nc,1);

% Now we can build the new coefficient matrix, and corresponding dpvar
Cnew = sparse(Cinew,Cjnew,Cval,nd*n,ntt*n);
ZQZ = dpvar(Cnew,degmat,varname,Q.dvarname,[n,n]);

end



function [ZQZ] = Z2QZ3construct(ZL,ZR,Q,n,var1,dum)
% This function computes the product bZLetas'*Q*bZRetath, producing a dpvar
% of dimensions [n,n];

nc = size(Q,2);
ZLdeg = ZL.degmat;
ZLntt = size(ZLdeg,1);
ZRdeg = ZR.degmat;
ZRntt = size(ZRdeg,1);

varname = [ZL.varname;dum.varname];

% We first build the new degmat
%row = repmat((1:ZRntt)',ZLntt,1);           % row indices of elements of Q
%col = kron((1:ZLntt)',ones(ZRntt,1));       % column indices

ntt = ZLntt * ZRntt;    nd = length(Q.dvarname)+1;
% degmat of the product Z2etas'*Q*Z2etath, accounting for the variable
% substitutions Z2sth --> Z2etas and Z2sth --> Z2etath
if strcmp(var1.varname{1},ZL.varname{1}) % var1 (s) is in the first position
    %degmat = [ZLdeg(row,2),zeros(length(row),1),ZLdeg(row,1)]  + ...
    %            [zeros(length(col),1),ZRdeg(col,2:-1:1)];
    degmat = [repmat(ZLdeg(:,2),ZRntt,1),zeros(ntt,1),repmat(ZLdeg(:,1),ZRntt,1)] + ...
                [zeros(ntt,1),kron(ZRdeg(:,2:-1:1),ones(ZLntt,1))];
else % var1 (s) is in the second position
    %degmat = [ZLdeg(row,1),zeros(length(row),1),ZLdeg(row,2)]  + ...
    %            [zeros(length(col),1),ZRdeg(col,:)];
    degmat = [repmat(ZLdeg(:,1),ZRntt,1),zeros(ntt,1),repmat(ZLdeg(:,2),ZRntt,1)] + ...
                [zeros(ntt,1),kron(ZRdeg,ones(ZLntt,1))];
end
%ntt = size(degmat,1);   nd = length(Q.dvarname)+1;

% Extract the coefficients and their locations from Q (Cval will be just
% ones)
[Ci,~,Cval] = find(Q.C);

% Each column in Q.C must be shifted to multiply with the appropriate row
% in the degmat
jndx_mat = reshape((1:ntt*n)',ZRntt,n*ZLntt);   % matrix of row indices in repeated degmat coefficients Q(1:ZLntt,1:end) should multiply with
jndx_mat = repmat(jndx_mat,[n,1]);              % matrix of indices indicating which row in degmat each coefficient Q(1:end,1:end) should be multiplied with
Cjnew = vec(jndx_mat);                          % adjusted column indices for coefficients in C

% Different coefficients Q(i,j) will add to the same element of the final
% object, we have to compress the row-indices accordingly
i_cmprss = repmat((0:nd:nd*(ZLntt-1))',n,1) + kron((0:nd*(ZLntt-1):(n-1)*nd*(ZLntt-1))',ones(ZLntt,1));
Cinew = Ci - repmat(i_cmprss,nc,1);

% Now we can build the new coefficient matrix, and corresponding dpvar
Cnew = sparse(Cinew,Cjnew,Cval,nd*n,ntt*n);
ZQZ = dpvar(Cnew,degmat,varname,Q.dvarname,[n,n]);

end


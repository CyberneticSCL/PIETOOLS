function [prog,Pop,LLL] = poslpivar(prog,n,I,d,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop,LLL] = poslpivar(prog,n,I,d,options) declares 
% the positive 4-PI operator. 
% P = T11 int(gs,s,a,b)
% Q(s) = g(s) T12 Z1(s) + int( gth T13 Z2(th,s) dth, s, b) + int( gth T14 Z3(th,s) dth, a, s)
% R0(s) = gs Z1(s) T22 Z1(s)
% R1(s,th) = gs Z1(s) T23 Z2(s,th) + gth Z3(th,s) T42 Z1(th) + 
%            int(geta Z2(eta,s)' T33 Z2(eta,th) deta,s,b)+int(geta Z3(eta,s) T43 Z2(eta,th) deta,th,s)
%           +int(geta Z3(eta,s) T44 Z3(eta,th) deta,a,th)
% R2(s,th) = R1(th,s)'
% T = [ T_{11}  T_{12} T_{13} T_{14}]
%     [ T_{21}  T_{22} T_{23} T_{24}] >0
%     [ T_{31}  T_{32} T_{33} T_{34}]
%     [ T_{41}  T_{42} T_{43} T_{44}]
% 
% where Z(x)= Z_d1(x) \otimes I_n and Z_d(x) is the vector of monomials in
% variables x of degree d1 or less. Z(x,y) = Z_{d1}(x) \otimes Z_{d2}(y)
% \otimes I_n. If the application is stability of time-delay systems, d1
% will eventually be more or less the degree of the multiplier and d2 more
% or less the degree of the kernel function.
% 
% INPUT 
%   prog: SOS program to modify.
%   n(1): dimension of real part
%   n(2): dimension of L2 part
%   I = [l u] interval of integration
%   -Optional INPUTS
%   d{1}: degree of s in Z1(s)
%   d{2}(1): degree of s in Z2(s,th), defaults to d(1)
%   d{2}(2): degree of th in Z2(s,th), defaults to d(1)
%   d{2}(3): joint degree of s,th in Z2(s,th), defaults to d(2,1)+d(2,2)
%   d{3}(1): degree of s in Z3(s,th), defaults to d(1)
%   d{3}(2): degree of th in Z3(s,th), defaults to d(1)
%   d{3}(3): joint degree of s,th in Z3(s,th), defaults to d(3,1)+d(3,2)
%   options.psatz=1 if this is a psatz term. options.psatz=0 otherwise
%   options.exclude is a length 4 binary vector where 
%      options.exclude(i)=1 if we want to set $T_{ij}=0$ for j=1...4
%   options.sep=1 if the kernel is separable, i.e. R1 = R2 
% 
% OUTPUT 
%   prog: modified SOS program with new variables and constraints
%   Pop: operator structure
%   - OPTIONAL OUTPUT
%   LLL: positive matrix variable used to create Pop
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - poslpivar
%
% Copyright (C)2019  M. Peet, S. Shivakumar
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS  - 7_26_2019
% MP 8_20_2019 - modified to add exclude option, input cell degrees for Z2,Z3 separately
% SS 8/26/2019 - added functionality to choose separable kernels, i.e.
% R1=R2
% SS 8/30/2019 - removed polynomial multiplications in Ri's to speed up the
% code

switch nargin
    case 2
        error(['Not enough inputs!'])
    case 3
        d = {1,[1,1,1],[1,1,1]};
        options.psatz=0;
        options.exclude=[0 0 0 0];
        options.diag=0;
        options.sep =0;
    case 4
        options.psatz=0;
        options.exclude=[0 0 0 0];
        options.diag=0;
        options.sep =0;
    case 5
        if ~isfield(options,'psatz')
            options.psatz=0;
        end
        if ~isfield(options,'exclude')
            options.exclude=[0 0 0 0];
        end
        if ~isfield(options,'diag')
            options.diag=0;
        end
        if ~isfield(options,'sep')
            options.sep=0;
        end
end
if length(n)~=2
    error('n must be a length 2 vector')
end

if I(1)>=I(2)
    error('I(1) must be less than I(2)')
end

if ~iscell(d)
    error('d is a 3-cell structure')
end

if length(d(:))==1
    d{2}=[d{1},d{1},2*d{1}];
    d{3}=d{2};
elseif length(d(:))==2
    if length(d{2})==1
        d{2}(2)=d{2}(1);
        d{2}(3)=d{2}(1);
    elseif length(d{2})==2
        d{2}(3)=d{2}(1);
    end
    d{3}=d{2};
else
    if length(d{2})==1
        d{2}(2)=d{2}(1);
        d{2}(3)=d{2}(1);
    elseif length(d{2})==2
        d{2}(3)=d{2}(1);
    end
    if length(d{3})==1
        d{3}(2)=d{3}(1);
        d{3}(3)=d{3}(1);
    elseif length(d{3})==2
        d{3}(3)=d{3}(1);
    end
end

if length(I)~=2
    error('I must be a length 2 vector')
end

pvar s theta;
var1 = s; var2 = theta;


n1=n(1);n2=n(2);
if n1==0
    if n2==0
        error(['Error in posopvar: All dimensions are zero'])
    end
end


% Sorting
% ExcludeL is a length-4 binary vector of terms to exclude
excludeL=options.exclude;

% in separable case R1=R2, so Ti3=Ti4 and T33=T44. We can reduce complexity
% by leveraging this relation. We will exclude Ti4 entirely. Use Ti3
% instead along with changing limits of integration. 8/26 - sachin
if options.sep ==1 
    excludeL(4)=1;
end

if n1==0
    excludeL(1)=1;
end
if n2==0
    excludeL([2:4])=[1 1 1];
end

% internal integration variable (replaces eta)
sss=polynomial(1,1,{'sss'},[1 1]);

% this defines the multiplier to be used later
if options.psatz==1
    gs=(var1-I(1))*(I(2)-var1);
    gth=(var2-I(1))*(I(2)-var2);
    geta=(sss-I(1))*(I(2)-sss);
else
    gs=polynomial(1);
    gth=polynomial(1);
    geta=polynomial(1);
end


if 0 
Z1s = Z{1};
Z2sth = Z{2};
Z3sth = Z{3};
nZ1=length(Z1s.degmat);
nZ2=length(Z2sth.degmat);
nZ3=length(Z3sth.degmat);
else
    % Constructing Z1(s)
nZ1=d{1}+1;
Z1degmat = [0:d{1}]';
Z1coeff = speye(nZ1);
Z1varname = var1.varname;
Z1matdim = [nZ1 1];
Z1s=polynomial(Z1coeff,Z1degmat,Z1varname,Z1matdim);

% Constructing Z2(s,th)
% In this implementation, Z2 will have degree d{2,2} in th and degree d{2,1} in s
% and max degree of s+th is d{2,3}. Similarly for Z3(s,th)



Z2degmat = [repmat([0:d{2}(1)]',d{2}(2)+1,1),vec(repmat(0:d{2}(2),d{2}(1)+1,1))];
Z2degmat(sum(Z2degmat,2)>d{2}(3),:)= [];
nZ2=size(Z2degmat,1);
Z2coeff = speye(nZ2);
Z2varname = [var1.varname; var2.varname];
Z2matdim = [nZ2 1];
Z2sth=polynomial(Z2coeff,Z2degmat,Z2varname,Z2matdim);

Z3degmat = [repmat([0:d{3}(1)]',d{3}(2)+1,1),vec(repmat(0:d{3}(2),d{3}(1)+1,1))];
Z3degmat(sum(Z3degmat,2)>d{3}(3),:)= [];
nZ3=size(Z3degmat,1);
Z3coeff = speye(nZ3);
Z3varname = [var1.varname; var2.varname];
Z3matdim = [nZ3 1];
Z3sth=polynomial(Z3coeff,Z3degmat,Z3varname,Z3matdim);
end


nBZ1=n2*nZ1;
nBZ2=n2*nZ2;
nBZ3=n2*nZ3;


bZ1s=[];
for i=1:n2
    bZ1s=blkdiag(bZ1s,Z1s);
end

bZ2sth=[];
for i=1:n2
    bZ2sth=blkdiag(bZ2sth,Z2sth);
end

bZ3sth=[];
for i=1:n2
    bZ3sth=blkdiag(bZ3sth,Z3sth);
end

bZ1th=subs(bZ1s,var1,var2);

bZ2ths=var_swap(bZ2sth,var1,var2);
bZ2etath=subs(bZ2sth,var1,sss);
bZ2etas=subs(bZ2ths,var2,sss);

bZ3ths=var_swap(bZ3sth,var1,var2);
bZ3etath=subs(bZ3sth,var1,sss);
bZ3etas=subs(bZ3ths,var2,sss);

% We now declare the positive matrix variable. 
% first compute the size of the matrix
dimL1=n1;
dimL2=(1-excludeL(2))*nBZ1;
dimL3=(1-excludeL(3))*nBZ2;
dimL4=(1-excludeL(4))*nBZ3;
dimLLL=dimL1+dimL2+dimL3+dimL4;
[prog,LLL]=sosposmatr(prog,dimLLL);


% Lets start by just assuming LLL has been partitioned
% In this case, things are more or less a repeat of the positive matrix
% variable case in sosposmatrvar.m with P replaced by LLL1 and Z with ZZZth
% Note, however that because this rearranges the order of the elements of
% ZZZth, we must alter all the manipulations or the result will be invalid.

ind{1}=1:n1; 
ind{2}=(dimL1+1):(dimL1+dimL2); 
ind{3}=(dimL1+dimL2+1):(dimL1+dimL2+dimL3);
ind{4}=(dimL1+dimL2+dimL3+1):dimLLL;


includeL = find(excludeL==0);


for i=includeL
    for j=includeL
        if i>=j         
            Q{i,j}=LLL(ind{i},ind{j});
            if i~=j
                Q{j,i}=Q{i,j}.';
            end
        end
    end
end

% This option is really only intended for time-delay systems, wherein we
% know these terms are zero and can be omitted
if options.diag==1 && options.exclude(2)~=1 % if the diagonal override is used, make $Q22$ block diagonal
    Q22temp=Q{2,2}; %this is of dimension nBZ1=n2*nZth=K*diag*nZth
    K=n2/diag; % this is the number of blocks
    for i=1:K
        irange=(diag*(i-1)*nZs+1):(i*diag*nZs);    % n2/
        Q22temp(irange,irange)=zeros(diag*nZs); % parts of Q22 which are not Zero
    end
    for i=1:K
        irange=(diag*(i-1)*nZs+1):(i*diag*nZs);    % n2/
        for j=1:K
            if i~=j
                jrange=(diag*(j-1)*nZs+1):(j*diag*nZs);    % n2/
                Q{2,2}(irange,jrange)=zeros(length(irange),length(jrange));
                % eliminate these terms so fewer variables appear in the polynomial manipulations
            end
        end
    end
    prog=sosmateq(prog,Q22temp); %constrain the off-diagonal terms to be zero
end

P=polynomial(zeros(n1));
QT=polynomial(zeros(n1,n2));
R0=polynomial(zeros(n2));
R1=R0;
R2=R0;


if excludeL(1)==0
    P=P+Q{1,1}*int(gs,var1,I(1),I(2));
    if excludeL(2)==0
        QT=QT+gs*Q{1,2}*bZ1s;
    end
    if excludeL(3)==0 && ~options.sep==1 % 8/26/19 changed to include R1=R2 - sachin
        QT=QT+int(gth*Q{1,3}*bZ2ths,var2,var1,I(2));
    elseif excludeL(3)==0 && options.sep==1
        QT=QT+int(gth*Q{1,3}*bZ2ths,var2,I(1),I(2));
    end
    if excludeL(4)==0
        QT=QT+int(gth*Q{1,4}*bZ2ths,var2,I(1),var1);
    end
end

if excludeL(2)==0
    R0 = constructR0(Q{2,2},Z1s,gs,n2,nZ1); 
    if excludeL(3)==0 && ~options.sep==1
        R1=R1+constructA1(Q{2,3},Z1s,Z2sth,gs,n2,nZ1,nZ2,var1,var2); 
    elseif excludeL(3)==0 && options.sep==1
        R1=R1+constructA1(Q{2,3},Z1s,Z2sth,gs,n2,nZ1,nZ2,var1,var2)+constructB1(Q{3,2},Z1s,Z2sth,gth,n2,nZ1,nZ2,var1,var2); % 8/30/2019
    end
    if excludeL(4)==0
        R1=R1+constructB1(Q{4,2},Z1s,Z2sth,gth,n2,nZ1,nZ2,var1,var2); 
    end
end

if excludeL(3)==0 && ~options.sep==1
        R1=R1+int(constructCA1(Q{3,3},Z2sth,geta,n2,nZ2,var1,var2,sss),sss,var1,I(2)); 
    if excludeL(4)==0
        R1=R1+int(constructCA1(Q{4,3},Z2sth,geta,n2,nZ2,var1,var2,sss),sss,var2,var1); 
    end
elseif excludeL(3)==0 && options.sep==1
    R1=R1+int(constructCA1(Q{3,3},Z2sth,geta,n2,nZ2,var1,var2,sss),sss,I(1),I(2)); 
end
if excludeL(4)==0
    R1=R1+int(constructCA1(Q{4,4},Z2sth,geta,n2,nZ2,var1,var2,sss),sss,I(1),var2);
end

if options.sep==1
    R2 = R1;
else
R2 = var_swap(R1,var1,var2).';
end

opvar Pop;
Pop.P = P;
Pop.Q1 = QT;
Pop.Q2 = QT'; Pop.R.R0 = R0; Pop.R.R1 = R1; Pop.R.R2 = R2;

Pop.dim = [n',n'];
Pop.I = I;
Pop.var1 = var1;
Pop.var2 = var2;
end

function [R0] = constructR0(Q11,ZZZth,g,n,nZth)
%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct M

avarname=[Q11.varname; ZZZth.varname];

bZdegmat=repmat(ZZZth.degmat,n,1);

[PIlist, PJlist]=find(Q11.coefficient);  
[PIlist2,ttt]=unique(PIlist);
PJlist2=PJlist(ttt,:);

[Prowu,Pcolu] = ind2sub([n*nZth n*nZth],PJlist2); % this returns the matrix locations for every term in P
adegmat = [Q11.degmat(PIlist2,:) bZdegmat(Prowu,:)+bZdegmat(Pcolu,:)];

[row,col] = ind2sub([n*nZth n*nZth],PJlist);
newrow=ceil(row/nZth);
newcol=ceil(col/nZth);
newidx=sub2ind([n n],newrow, newcol);

coeff=sparse(PIlist,newidx,1);
amatdim=[n n];
R0=g*polynomial(coeff,adegmat,avarname,amatdim);   
end

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

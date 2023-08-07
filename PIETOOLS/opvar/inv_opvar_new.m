function [Pinv] = inv_opvar_new(Pop, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pinv] = inv_opvar(Pop, tol) computes the inverse operator of the
% operator
% INPUT
%   Pop: positive definite opvar to invert
%   tol: order of coefficients in the polynomials that should be truncated
%
% OUTPUT
%   Pinv: inverse opvar object. Inverse opvar is a numerical inversion and
%   should be used with care and reservations.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - inv_opvar
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
% Initial coding SS - 5_31_2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if nargin==1
    tol=1e-7;
end

if ~isa(Pop,'opvar')|| any(Pop.dim(:,1)~=Pop.dim(:,2))
    error('Only opvar class objects with equal input and output dimensions can be inverted using this function');
end

Pop = clean(Pop,tol);

if all(isequal(Pop.R.R1,Pop.R.R2)) % if pi operator is separable then use analytical formulae from the old code
    Pinv = inv_opvar_old(Pop);
    return;
end


X = Pop.I; var1 = Pop.var1; var2 = Pop.var2;
opvar Pinv;
Pinv.I = X; Pinv.var1 = var1; Pinv.var2 = var2;

if all(Pop.dim(2,:)==0)
    if isa(Pop.P,'polynomial')
        Pop.P = double(Pop.P);
    end
    Pinv.P = inv(Pop.P);
elseif all(Pop.dim(1,:)==0)
    %     [F,G,R0inv] = getsemisepmonomials(Pop);
    %     F{1} = - F{1}; F{2} = - F{2};
    %
    %     A = subs([G{1}*F{1} G{1}*F{2}; -G{2}*F{1} -G{2}*F{2}],var2,var1);
    %
    %     U = polynomial(eye(size(A)));
    %     Uk = U;

    % find max terms needed for approximation
    %     normA = sqrt(max(abs(eig(double(int(A'*A,var1,X(1),X(2)))))));

    orderapp = 5;
    [Fd1, Fd2, Gd1, Gd2, Ad, R0inv, Rinvd, dx] = getsemisepmonomials_discrete(Pop);
     
    for i=1:size(Ad,1)
        tmp(i,:,:) = squeeze(Ad(i,:,:))'*squeeze(Ad(i,:,:));
        Ud(i,:,:) = eye(size(Ad,2),size(Ad,3));
    end
    Udinv = Ud;
    Udk = Ud; Udinvk = Udinv;


    normA = sqrt(max(abs(eig(trapz(dx,tmp)))));
    Nmax = 1;
    while (normA^Nmax)/factorial(Nmax)>tol
        Nmax = Nmax+1;
    end


    for i=1:Nmax
        % remove terms from integrand that would be too small after Nmax
        % integrations
        for i=1:size(Ad,1)
            tmp(i,:,:) = squeeze(Ad(i,:,:))*squeeze(Udk(i,:,:));
        end
        Udk = trapz(dx,tmp);
        Ud = Ud+Udk;
    end

    tmp = squeeze(U(end,:,:));

    U11 = tmp(1:size(G1,2), 1:size(F1,3));
    U12 = tmp(1:size(G1,2), size(F1,3)+1:end);
    U21 = tmp(size(G1,2)+1:end, 1:size(F1,3));
    U22 = tmp(size(G1,2)+1:end, size(F1,3)+1:end);


    if rcond(U22)<eps
        error('Given PI operator is likely non-invertible');
    else
        P = [zeros(size(U11)) zeros(size(U12)); U22\U21 eye(size(U22))];

        for i=1:size(Ad,1)
            C(i,:,:) = [Fd1(i,:,:) Fd2(i,:,:)]; B(i,:,:) = [Gd1(i,:,:); -Gd2(i,:,:)];
        end
        
       

        % finding U-inverse
        for i=1:Nmax
            % remove terms from integrand that would be too small after Nmax
            % integrations
            for i=1:size(Ad,1)
                tmp(i,:,:) = squeeze(Udink(i,:,:))*squeeze(Ad(i,:,:));
            end
    
            Udinvk = trapz(dx,tmp);
            Udinv = Udinv-Udinvk;
        end

        Pinv.R.R0 = R0inv;

        for i=1:size(Ad,1)
                tmp(i,:,:) = squeeze(C(i,:,:))*squeeze(U(i,:,:))*(eye(size(P))-P)*squeeze(Uinv(i,:,:))*squeeze(B(i,:,:));
        end
        tmp2 = polynomial(zeros(size(tmp,[2,3])));
        for i=1:size(tmp,2)
            for j=1:size(tmp,3)
                Data1=squeeze(tmp(:,i,j))';
                tempCoeffs =polyfit([X(1):dx:X(2)],Data1,orderapp); % uses matlab internal polynomial representation
                tmp2(i,j)=var1.^(orderapp:-1:0)*tempCoeffs';
            end
        end

        Pinv.R.R1 = tmp2*subs(R0inv,var1,var2);
        tmp = -C*U*P*Uinv*B*subs(R0inv,var1,var2);
        Pinv.R.R2 = tmp;
        Pinv = clean(Pinv,tol);
    end
else
    opvar A B C D;
    Ad.I = X; Ad.var1 = Pop.var1; Ad.var2 = Pop.var2;
    B.I = X; B.var1 = Pop.var1; B.var2 = Pop.var2;
    C.I = X; C.var1 = Pop.var1; C.var2 = Pop.var2;
    D.I = X; D.var1 = Pop.var1; D.var2 = Pop.var2;
    Ad.P = Pop.P; B.Q1 = Pop.Q1; C.Q2 = Pop.Q2; D.R = Pop.R;
    %     Dinv = opvar_inverse_iterative(D,tol);
    Ainv = inv_opvar(Ad,tol);
    %     TA = opvar_inverse_iterative(A-B*Dinv*C,tol);
    TB = inv_opvar(D-C*Ainv*B,tol);
    %     Pinv = [TA,         -Ainv*B*TB;
    %             -Dinv*C*TA, TB];

    Pinv = [Ainv+Ainv*B*TB*C*Ainv, -Ainv*B*TB; -TB*C*Ainv, TB];
    Pinv = clean(Pinv,tol);
end
end
function [F1,F2,G1,G2,A,Rinv, Rtemp,dX] = getsemisepmonomials_discrete(P)
if ~isa(P,'opvar')
    error('Input should be an opvar class object');
end
X = P.I; var1 = P.var1; var2 = P.var2;

% first separate polynomials in s and theta
Ra = P.R.R1;
Rb = P.R.R2;

% Error check: Change Rinv, Ra, Rb to polynomials if they are not polynomials
R1temp = polynomial(Ra);
R2temp = polynomial(Rb);

% Error check: fix degmats if var1, var2 are missing
if isempty(R1temp.degmat)% if neither s nor theta is present
    R1temp = polynomial(R1temp.coefficient,zeros(size(R1temp.degmat,1),2),[var1.varname;var2.varname],R1temp.matdim);
elseif size(R1temp.degmat,2)<2 % if only one of s or theta is present
    if ismember(R1temp.varname,var1.varname)
        missingVar = var2.varname;
    else
        missingVar = var1.varname;
    end
    R1temp = polynomial(R1temp.coefficient,[R1temp.degmat,zeros(size(R1temp.degmat,1),1)],[R1temp.varname;missingVar],R1temp.matdim);
end
if isempty(R2temp.degmat) % if neither s nor theta is present
    R2temp = polynomial(R2temp.coefficient,zeros(size(R2temp.degmat,1),2),[var1.varname;var2.varname],R2temp.matdim);
elseif size(R2temp.degmat,2)<2 % if only s or theta is present
    if ismember(R2temp.varname,var1.varname)
        missingVar = var2.varname;
    else
        missingVar = var1.varname;
    end
    R2temp = polynomial(R2temp.coefficient,[R2temp.degmat,zeros(size(R2temp.degmat,1),1)],[R2temp.varname;missingVar],R2temp.matdim);
end
Ra = R1temp; Rb = R2temp;

Ra_vnames = Ra.varname;
Rb_vnames = Rb.varname;

Ra_maxdeg = max(Ra.degmat,[],1);
Rb_maxdeg = max(Rb.degmat,[],1);

var1_loc_Ra = ismember(Ra_vnames,var1.varname);
var2_loc_Ra = ismember(Ra_vnames,var2.varname);
if Ra_maxdeg(var1_loc_Ra)<= Ra_maxdeg(var2_loc_Ra)
    [newdegmat,idx] = sortrows(Ra.degmat);
    [val,~,~] = unique(full(Ra.degmat(:,var1_loc_Ra)));
    Ra = polynomial(Ra.coefficient(idx,:),Ra.degmat(idx,:),Ra.varname,Ra.matdim);
    F{1} = kron(var1.^(val'),eye(size(Ra)));
    G{1} = {};
    for i=1:length(val)
        loctemp = find(newdegmat(:,var1_loc_Ra)==val(i));
        temp = polynomial(Ra.coefficient(loctemp,:),Ra.degmat(loctemp,var2_loc_Ra),var2.varname,Ra.matdim);
        G{1} = [G{1};temp];
    end
else
    [newdegmat,idx] = sortrows(Ra.degmat,2);
    [val,~,~] = unique(full(Ra.degmat(:,var2_loc_Ra)));
    Ra = polynomial(Ra.coefficient(idx,:),Ra.degmat(idx,:),Ra.varname,Ra.matdim);
    G{1} = kron(var2.^val,eye(size(Ra)));
    F{1} = {};
    for i=1:length(val)
        loctemp = find(newdegmat(:,var2_loc_Ra)==val(i));
        temp = polynomial(Ra.coefficient(loctemp,:),Ra.degmat(loctemp,var1_loc_Ra),var1.varname,Ra.matdim);
        F{1} = [F{1},temp];
    end
end

var1_loc_Rb = ismember(Rb_vnames,var1.varname);
var2_loc_Rb = ismember(Rb_vnames,var2.varname);
if Rb_maxdeg(var1_loc_Rb)<= Rb_maxdeg(var2_loc_Rb)
    [newdegmat,idx] = sortrows(Rb.degmat);
    [val,~,~] = unique(full(Rb.degmat(:,var1_loc_Ra)));
    Rb = polynomial(Rb.coefficient(idx,:),Rb.degmat(idx,:),Rb.varname,Rb.matdim);
    F{2} = kron(var1.^(val'),eye(size(Rb)));
    G{2} = {};
    for i=1:length(val)
        loctemp = find(newdegmat(:,var1_loc_Rb)==val(i));
        temp = polynomial(Rb.coefficient(loctemp,:),Rb.degmat(loctemp,var2_loc_Rb),var2.varname,Rb.matdim);
        G{2} = [G{2};temp];
    end
else
    [newdegmat,idx] = sortrows(Rb.degmat,2);
    [val,~,~] = unique(full(Rb.degmat(:,var2_loc_Rb)));
    Rb = polynomial(Rb.coefficient(idx,:),Rb.degmat(idx,:),Rb.varname,Rb.matdim);
    G{2} = kron(var2.^val,eye(size(Rb)));
    F{2} = {};
    for i=1:length(val)
        loctemp = find(newdegmat(:,var2_loc_Rb)==val(i));
        temp = polynomial(Rb.coefficient(loctemp,:),Rb.degmat(loctemp,var1_loc_Rb),var1.varname,Rb.matdim);
        F{2} = [F{2},temp];
    end
end

% finding U-inverse
N=100; orderapp=max(P.R.R0.degmat(:));
dX = linspace(X(1),X(2),N);
dx = dX(2)-dX(1);
ii=0;
Rtemp = zeros(N,size(P.R.R0,1),size(P.R.R0,2));
for ss=dX
    ii=ii+1;
    Rtemp(ii,:,:)= inv(double(subs(P.R.R0,var1,ss))); % Calculates the value of the inverse of S at every point in the interval
    tmp1 = squeeze(Rtemp(ii,:,:))*subs(F{1},var1,ss);
    tmp2 = squeeze(Rtemp(ii,:,:))*subs(F{2},var1,ss); 
    tmp3 = subs(G{1},var2,ss); tmp4 = subs(G{2},var2,ss);
    
    F1(ii,:,:)= tmp1;  F2(ii,:,:)= tmp2; G1(ii,:,:)= tmp3; G2(ii,:,:)= tmp4;
    A(ii,:,:)= [tmp3*tmp1, tmp3*tmp2; -tmp4*tmp1, -tmp4*tmp2];
end
Rinv = polynomial(zeros(size(P.R.R0)));
for i=1:size(P.R.R0,1)
    for j=1:size(P.R.R0,2)
        Data1=squeeze(Rtemp(i,j,:))';
        tempCoeffs =polyfit([X(1):dx:X(2)],Data1,orderapp); % uses matlab internal polynomial representation
        Rinv(i,j)=var1.^(orderapp:-1:0)*tempCoeffs';
    end
end

end

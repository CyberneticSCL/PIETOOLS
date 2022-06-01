function [Pinv] = inv_opvar(Pop, tol)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    tol=1e-8;
end

if ~isa(Pop,'opvar')|| ~any(Pop.dim(:,1)==Pop.dim(:,2))
    error('Only symmetrix opvar class objects can be inverted using this function');
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
    [F,G,R0inv] = getsemisepmonomials(Pop);
    F{1} = - F{1}; F{2} = - F{2};
    
    A = subs([G{1}*F{1} G{1}*F{2}; -G{2}*F{1} -G{2}*F{2}],var2,var1);
    
    U = polynomial(eye(size(A)));
    Uk = U;
    
    % find max terms needed for approximation
    normA = sqrt(max(abs(eig(double(int(A'*A,var1,X(1),X(2)))))));
    Nmax = 1;
    while (normA^Nmax)/factorial(Nmax)>eps
        Nmax = Nmax+1;
    end
    for i=1:Nmax
        Uk = int(A*Uk,var1,X(1),var1);
        if isa(Uk, 'double')
            Uk = polynomial(Uk);
        end
        if max(abs(Uk.coefficient(:)))<tol
            break;
        end
        U = U+Uk;
    end
    
    U.coefficient(find(abs(U.coefficient)<tol)) = 0;
    
    
    U11 = U(1:size(G{1},1), 1:size(F{1},2));
    U12 = U(1:size(G{1},1), size(F{1},2)+1:end);
    U21 = U(size(G{1},1)+1:end, 1:size(F{1},2));
    U22 = U(size(G{1},1)+1:end, size(F{1},2)+1:end);
    U21 = double(subs(U21,var1,X(2)));
    U22 = double(subs(U22,var1,X(2)));
    
    
    
    if rcond(U22)<eps
        error('Given PI operator is likely non-invertible');
    else
        C = [F{1} F{2}]; B = [G{1}; -G{2}];
        P = [zeros(size(U11)) zeros(size(U12)); U22\U21 eye(size(U22))];
        
        % finding U-inverse
        Uinv = polynomial(eye(size(A)));
        Ukinv = Uinv;
        for i=1:Nmax
            Ukinv = int(Ukinv*A,var1,X(1),var1);
            if isa(Ukinv, 'double')
                Ukinv = polynomial(Ukinv);
            end
            if max(abs(Ukinv.coefficient(:)))<tol
                break;
            end
            Uinv = Uinv-Ukinv;
        end
        
        Uinv.coefficient(find(abs(Uinv.coefficient)<tol)) = 0;
        Uinv = subs(Uinv,var1,var2);
        
        Pinv.R.R0 = R0inv;
        Pinv.R.R1 = C*U*(eye(size(P))-P)*Uinv*B*subs(R0inv,var1,var2);
        Pinv.R.R2 = -C*U*P*Uinv*B*subs(R0inv,var1,var2);
    end
else
    opvar A B C D;
    A.I = X; A.var1 = Pop.var1; A.var2 = Pop.var2;
    B.I = X; B.var1 = Pop.var1; B.var2 = Pop.var2;
    C.I = X; C.var1 = Pop.var1; C.var2 = Pop.var2;
    D.I = X; D.var1 = Pop.var1; D.var2 = Pop.var2;
    A.P = Pop.P; B.Q1 = Pop.Q1; C.Q2 = Pop.Q2; D.R = Pop.R;
%     Dinv = opvar_inverse_iterative(D,tol);
    Ainv = inv_opvar(A,tol);
%     TA = opvar_inverse_iterative(A-B*Dinv*C,tol);
    TB = inv_opvar(D-C*Ainv*B,tol);
%     Pinv = [TA,         -Ainv*B*TB;
%             -Dinv*C*TA, TB];
        
    Pinv = [Ainv+Ainv*B*TB*C*Ainv, -Ainv*B*TB; -TB*C*Ainv, TB];
end
end
function [F,G, Rinv] = getsemisepmonomials(P)
if ~isa(P,'opvar')
    error('Input should be an opvar class object');
end
X = P.I; var1 = P.var1; var2 = P.var2;

% finding U-inverse  
N=100; orderapp=4; 
dx = (X(2)-X(1))/N; 
ii=0;
Rtemp = zeros(size(P.R.R0,1),size(P.R.R0,2),N);
for ss=[X(1):dx:X(2)]
    ii=ii+1;
    Rtemp(:,:,ii)= inv(double(subs(P.R.R0,var1,ss))); % Calculates the value of the inverse of S at every point in the interval
end
Rinv = polynomial(zeros(size(P.R.R0)));
for i=1:size(P.R.R0,1)
    for j=1:size(P.R.R0,2)
        Data1=squeeze(Rtemp(i,j,:))';
        tempCoeffs =polyfit([X(1):dx:X(2)],Data1,orderapp); % uses matlab internal polynomial representation 
        Rinv(i,j)=var1.^(orderapp:-1:0)*tempCoeffs';
    end
end

Rinv.coefficient(find(abs(Rinv.coefficient)<1e-10)) = 0;
Ra = Rinv*P.R.R1;
Rb = Rinv*P.R.R2;

% Error check: Change Rinv, Ra, Rb to polynomials if they are not polynomials
R0temp = polynomial(Rinv);
R1temp = polynomial(Ra);
R2temp = polynomial(Rb);

% Error check: fix degmats if var1, var2 are missing
if isempty(R0temp.degmat)% if s is not present
    R0temp = polynomial(R0temp.coefficient,zeros(size(R0temp.degmat,1),1),var1.varname,R0temp.matdim);
end
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
Rinv = R0temp; Ra = R1temp; Rb = R2temp;

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
end

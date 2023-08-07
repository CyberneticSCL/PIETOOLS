function [Pinv] = inv_opvar_new(Pop, tol,N)
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
    si = linspace(X(1),X(2),N); ds = si(2)-si(1);
    [F,G,R0inv] = getsemisepmonomials(Pop,N);
    F{1} = - F{1}; F{2} = - F{2};
    
    % R_1a = F1, R_1b = G1, R_2a = G2, R_2b = G2, R_0 = R0

    % find A(s_i) = [R_{1,b}(s_i)R_0(s_i)^{-1}R_{1,a}(s_i)&
    %                       R_{1,b}(s_i)R_0(s_i)^{-1}R_{2,a}(s_i)\\
    %                         -R_{2,b}(s_i)R_0(s_i)^{-1}R_{1,a}(s_i)&
    %                            -R_{2,b}(s_i)R_0(s_i)^{-1}R_{2,a}(s_i)]
    
    for i=1:length(si)
        sval = si(i);
        tmp1 = subs(G{1},var2,sval)*R0inv{i};
        tmp2 = -subs(G{2},var2,sval)*R0inv{i};
        A{i} = double([tmp1*subs(F{1},var1,sval),tmp1*subs(F{2},var1,sval);
                       tmp2*subs(F{1},var1,sval),tmp2*subs(F{2},var1,sval)]);
    end
    
    [U,V] = solve_volterra(A,si);
    
    tmp = U{end};

    U11 = tmp(1:size(G{1},1), 1:size(F{1},2));
    U12 = tmp(1:size(G{1},1), size(F{1},2)+1:end);
    U21 = tmp(size(G{1},1)+1:end, 1:size(F{1},2));
    U22 = tmp(size(G{1},1)+1:end, size(F{1},2)+1:end);

    if rcond(U22)<eps
        error('Given PI operator is likely non-invertible');
    else
        C = [F{1} F{2}]; B = [G{1}; -G{2}];
        P = [zeros(size(U11)) zeros(size(U12)); U22\U21 eye(size(U22))];
        
        for i=1:length(si)
            for j=1:length(si)
                Ct = double(subs(C,var1,si(i))); Bt = double(subs(B,var2,si(j)));
                R1{i,j}=R0inv{i}*Ct*U{i}*(eye(size(P))-P)*V{j}*Bt*R0inv{j};
                R2{i,j}=-R0inv{i}*Ct*U{i}*P*V{j}*Bt*R0inv{j};
            end
        end
        
        % fit polynomial
        orderapp=4;
        Rinv = polynomial(zeros(size(Pop.R.R0))); R1inv = Rinv; R2inv = Rinv;
        Rtemp = cat(3,R0inv{:}); 
        R1temp = reshape( cat(3,R1{:})  , [size(Pop.R.R0,2),size(Pop.R.R0,2),length(si),length(si)]); 
        R2temp = reshape( cat(3,R2{:})  , [size(Pop.R.R0,2),size(Pop.R.R0,2),length(si),length(si)]); 
%         for i=1:size(Pop.R.R0,1)
%             for j=1:size(Pop.R.R0,2)
%                 Data1=squeeze(Rtemp(i,j,:))';
%                 tempCoeffs =polyfit([X(1):dx:X(2)],Data1,orderapp); % uses matlab internal polynomial representation
%                 Rinv(i,j)=var1.^(orderapp:-1:0)*tempCoeffs';
%             end
%         end
        Pinv = {};
        Pinv{1} = R0inv;
        Pinv{2} = R1;
        Pinv{3} = R2;
%         Pinv = clean(Pinv,tol);
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
    Pinv = clean(Pinv,tol);
end
end
function [F,G,Rtemp] = getsemisepmonomials(P,N)
if ~isa(P,'opvar')
    error('Input should be an opvar class object');
end
X = P.I; var1 = P.var1; var2 = P.var2;


    % finding R0-inverse
    orderapp=max(P.R.R0.degmat(:));
    dx = linspace(X(1),X(2),N);
    for ss=1:length(dx)
        if poly2double(P.R.R0)
            Rtemp{ss} = inv(double(P.R.R0));
        else
            Rtemp{ss}= inv(double(subs(P.R.R0,var1,dx(ss)))); % Calculates the value of the inverse of S at every point in the interval
        end
end


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
end

function [U,V] = solve_volterra(A,si)
ds = si(end)-si(end-1);
U{1} = eye(size(A{1}));
V{1} = eye(size(A{1}));
bu= U{1}; bv = V{1};
for i=2:length(si)
    bu = bu + A{i-1}*U{i-1}*ds/2;
    Au = U{1}-A{i}*ds/2; 
    U{i} = Au\bu;
    bv = bv+V{i-1}*A{i-1}*ds/2;
%     V{i} = bv/Au;
    V{i} = inv(U{i});
end
end

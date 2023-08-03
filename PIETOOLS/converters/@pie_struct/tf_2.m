function varargout = tf_2(PIE,w)
% assuming inputs and outputs are finite dimensional
% pde doesnt have ode coupling

nx = PIE.T.dim(:,2);
nz = PIE.C1.dim(1,1);
nw = PIE.B1.dim(1,2);

C = PIE.C1; B = PIE.B1; Tw = PIE.Tw;
T = PIE.T; A = PIE.A;

opvar Z Zb;
Z.dim = PIE.C1.dim;
Zb.dim = PIE.B1.dim;

bigC = [C,Z;Z,C];
bigT = [-A,w*T; -w*T,-A];
bigB = [B-w*Tw,Zb;Zb,B-w*Tw];

Cq1 = bigC.Q1;
Bq2 = bigB.Q2;

[F,G] = getsemisepmonomials(bigT);
F{1} = - F{1}; F{2} = - F{2};

% R_1a = F1, R_1b = G1, R_2a = G2, R_2b = G2, 

% find A(s_i) = [R_{1,b}(s_i)R_0(s_i)^{-1}R_{1,a}(s_i)&
%                       R_{1,b}(s_i)R_0(s_i)^{-1}R_{2,a}(s_i)\\
%                         -R_{2,b}(s_i)R_0(s_i)^{-1}R_{1,a}(s_i)&
%                            -R_{2,b}(s_i)R_0(s_i)^{-1}R_{2,a}(s_i)]
N = 100;
si = linspace(X(1),X(2),N); ds = s(2)-s(1);
for i=1:length(si)
    sval = si(i);
    tmp = subs(R0,var1,sval);
    R0inv(i,:,:)=inv(tmp);
    tmp1 = subs(G{1},var2,sval)/squeeze(tmp);
    tmp2 = -subs(G{2},var2,sval)/squeeze(tmp);
    A(i,:,:) = [tmp1*subs(F{1},var1,sval),tmp1*subs(F{2},var1,sval);
                  tmp2*subs(F{1},var1,sval),tmp2*subs(F{2},var1,sval)];
    U(i,:,:) = eye(size(A(i,:,:),[2,3]));
end
Uk = U;
V = U; Vk = U;

for j=1:Nmax
for i=1:size(A,1)
    tmp(i,:,:) = squeeze(A(i,:,:))*squeeze(Uk(i,:,:));
    tmp2(i,:,:) = squeeze(Vk(i,:,:))*squeeze(A(i,:,:));
    Uk(i,:,:) = trapz(si(1:i),tmp(1:i,:,:));
    Vk(i,:,:) = trapz(si(1:i),tmp2(1:i,:,:));
end
U=U+Uk; V=V+Vk;
end

q = size(G{2},1); p = size(F{1},2);
tmp = squeeze(U(end,:,:));
U22 = tmp(p+1:end,p+1:end); U21 = tmp(p+1:end,1:p);
P = [zeros(p), zeros(p,q); U22\U21, eye(q)];
I_P = eye(size(P))-P;


for i=1:size(A,1)
    lC1(i,:,:) = subs(Cq1*squeeze(R0inv(i,:,:))*[F{1},F{2}]*squeeze(U(i,:,:))*IP,var1,si(i));    
    lC2(i,:,:) = subs(Cq1*squeeze(R0inv(i,:,:))*[F{1},F{2}]*squeeze(U(i,:,:))*P,var1,si(i));    
    rB(i,:,:) = subs(squeeze(V(i,:,:))*[G{1};-G{2}]*squeeze(R0inv(i,:,:))*Bq2,var1,si(i));
end
for i=1:size(A,1)
    TF(i,:,:) = subs(Cq1*squeeze(R0inv(i,:,:))*Bq2,var1,si(i));
    TF1(i,:,:) = trapz(si(i:end),squeeze(lC1(i:end,:,:))*squeeze(rB(i:end,:,:)));
    TF2(i,:,:) = trapz(si(1:i),squeeze(lC2(1:i,:,:))*squeeze(rB(1:i,:,:)));
end
TF = TF+TF1+TF2;
tmp = squeeze(trapz(si, TF)); % TF is size n_grid*2*nz*2*nw;

varargout{1} = tmp(1:nz,1:nw);
varargout{2} = tmp(1:nz,nw+1:2*nw);
end
function [F,G] = getsemisepmonomials(P)
if ~isa(P,'opvar')
    error('Input should be an opvar class object');
end
X = P.I; var1 = P.var1; var2 = P.var2;

if poly2double(P.R.R0)
    Rinv = inv(double(P.R.R0));
else
    % finding R0-inverse
    N=100; orderapp=max(P.R.R0.degmat(:));
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
end


Ra = P.R.R1;
Rb = P.R.R2;


R1temp = polynomial(Ra);
R2temp = polynomial(Rb);

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
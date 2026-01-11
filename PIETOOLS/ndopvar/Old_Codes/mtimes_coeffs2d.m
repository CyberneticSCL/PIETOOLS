function P_str = mtimes_coeffs2d(C_str,D_str)
% P_STR = MTIMES_COEFFS(C_STR,D_STR) computes the composition of the
% 3-PI operators defined by structs C_STR and D_STR, such that
%   P_STR = opvar2d2coeffs(coeffs2opvar2d(C_STR)*coeffs2opvar2d(D_STR));
%
% INPUTS
% - C_str,D_str:    structs representing opvar2d objects in the format
%                   specified in opvar2d2coeffs;
%
% OUTPUTS
% - P_str


% Check that coefficients are appropriately specified
dim1 = C_str.dim;   dim2 = D_str.dim;
if dim1(2)~=dim2(1)
    error("Inner dimensions of operators to compose should match.")
end
m = dim1(1);    q = dim1(2);    n = dim2(2);
d1 = C_str.deg;     d2 = D_str.deg;
if d1~=d2
    error("Degrees of the opvar2d objects must match.")
end
d = d1;
% C0 = C_str.C{1};    C1 = C_str.C{2};      C2 = C_str.C{3};
% D0 = D_str.C{1};    D1 = D_str.C{2};      D2 = D_str.C{3};
% if any(size(C1)~=size(C2)) || any(size(D1)~=size(D2))
%     error("Coefficients representing each parameter should be of same size.")
% end
% d = size(C1,1)/m-1;
% if size(C1,2)~=q*(d+1) || size(D1,1)~=q*(d+1) || size(D1,2)~=n*(d+1)
%     error("Coefficients representing the operators should be specified as m*(d+1) x n*(d+1) arrays.")
% end

% Check that the domain is appropriately specified
dom1 = C_str.dom;       dom2 = D_str.dom;
if ~isa(dom1,'double') || ~numel(dom1)==2
    error("Spatial domain should be specified as 1x2 array")
end
if any(dom1~=dom2)
    error("Spatial domains of the operators should match.")
end
a = dom1(1,1);       b = dom1(1,2);


% Declare the matrices (Iq o A) and (Iq o B), where
% B-A = int_{a}^{b} Zd(s)*Zd(s)^T ds 
d_arr = (0:d)+(0:d)'+1;
vals = 1./d_arr(:);
r_idcs = repmat((1:d+1)',[1,d+1]);
c_idcs = repmat((1:d+1),[d+1,1]);
%nd = numel(d_arr);
%r_idcs = r_idcs(:) + repmat((0:q-1)*(d+1),[nd,1]);
%c_idcs = c_idcs(:) + repmat((0:q-1)*(d+1),[nd,1]);
if a==0
    Amat = sparse(q*(d+1),q*(d+1));
elseif a==1
    Amat = spIkron(q,r_idcs,c_idcs,vals,d+1,d+1);
    %avals = repmat(vals,[q,1]);
    %Amat = sparse(r_idcs,c_idcs,avals,q*(d+1),q*(d+1));
else
    avals = vals.*(a.^d_arr(:));
    Amat = spIkron(q,r_idcs,c_idcs,avals,d+1,d+1);
    %avals = repmat(vals.*(a.^d_arr(:)),[q,1]);
    %Amat = sparse(r_idcs,c_idcs,avals,q*(d+1),q*(d+1));
end
if b==0
    Bmat = sparse(q*(d+1),q*(d+1));
elseif b==1
    Bmat = spIkron(q,r_idcs,c_idcs,vals,d+1,d+1);
    %bvals = repmat(vals,[q,1]);
    %Bmat = sparse(r_idcs,c_idcs,bvals,q*(d+1),q*(d+1));
else
    bvals = vals.*(b.^d_arr(:));
    Bmat = spIkron(q,r_idcs,c_idcs,bvals,d+1,d+1);
    %bvals = repmat(vals.*(b.^d_arr(:)),[q,1]);
    %Bmat = sparse(r_idcs,c_idcs,bvals,q*(d+1),q*(d+1));
end

% Declare the permutation matrix (Iq o Sd) for Sd such that
%   int Z_{d}(s) Z_{d}(s)' = Sd kron(eye(d+1),Z_{2d+1}(s))
r_idcs = repmat((1:d+1)',[1,d+1]);
c_idcs = (1:d+1)'+(1:d+1)+(2*d+2)*(0:d);
Sd = spIkron(q,r_idcs,c_idcs,vals,d+1,(d+1)*(2*d+2));
Sdt = spIkron(q,c_idcs,r_idcs,vals,(d+1)*(2*d+2),d+1);

% Declare the permuation matrices (In o Ed) and (Im o Ed)^T for Ed s.t.
%   kron(Z_{d}(s),Z_{2d+1}(s)) = Ed*Z_{3d+1}(s)
vals = ones((d+1)*(2*d+2),1);
r_idcs = (1:(d+1)*(2*d+2))';
c_idcs = (1:2*d+2)'+(0:d);
Ed = spIkron(n,r_idcs(:),c_idcs(:),vals,(d+1)*(2*d+2),(3*d+2));
Edt = spIkron(m,c_idcs(:),r_idcs(:),vals,(3*d+2),(d+1)*(2*d+2));

% Declare the permutation matrix (Im o Fd)^T for Fd such that
%   kron(Z_{d}(s),Z_{d}(s)) = Fd*Z_{3d+1}(s)
r_idcs = (1:d+1)' + (0:d);
c_idcs = (1:(d+1)^2);
Fdt = spIkron(m,r_idcs(:),c_idcs(:),1,3*d+2,(d+1)^2);

% Declare the permutation matrix (Iq o Hd) for Hd such that
%   Z_{d}(s)*Zd_{s}' = Hd*(I_{d+1} o Z_{3d+1}(s))
Hd_vals = ones((d+1)^2,1);
r_idcs = repmat((1:d+1)',[1,d+1]);
c_idcs = (1:d+1)' + (0:d)*(3*d+3);
Hd = spIkron(q,r_idcs(:),c_idcs(:),Hd_vals,(d+1),(d+1)*(3*d+2));

% Declare the permutation matrices (In o E1) and (Im o E1)^T such that
%    Z_{d}(s) = E1*Z_{3d+1}(s)
E1 = spIkron(n,(1:d+1),(1:d+1),1,d+1,3*d+2);
E1t = spIkron(m,(1:d+1),(1:d+1),1,3*d+2,d+1);


% Decompose the 2D operator into 3 1D operators
[C0,C1,C2] = get1Dop(C_str);
[D0,D1,D2] = get1Dop(D_str);
C0_I = kronI_coeffs(C0,d+1);
D0_I = kronI_coeffs(D0,3*d+2);
C12 = minus_coeffs(C1,C2);
D21 = minus_coeffs(D2,D1);
C12_I = kronI_coeffs(C12,2*d+2);
D21_I = kronI_coeffs(D21,2*d+2);




% Finally, compute the coefficients defining the composition of the
% operators associated with C and D

fctrC0 = mat_times_coeffs(Fdt,C0_I);
fctrD0 = mat_times_coeffs(Hd,D0_I);
CBD = mtimes_coeffs(C2,mat_times_coeffs(Bmat,D1));
CAD = mtimes_coeffs(C1,mat_times_coeffs(Amat,D2));
AB_diff = minus_coeffs(CBD,CAD);
trm1 = mat_times_coeffs(E1t,coeffs_times_mat(AB_diff,E1));
fctrD = mat_times_coeffs(Sd,coeffs_times_mat(D21_I,Ed));
fctrC = coeffs_times_mat(mat_times_coeffs(Edt,C12_I),Sdt);

fctrC_new = plus_coeffs(fctrC0,fctrC);
fctrD_new = plus_coeffs(fctrD0,fctrD);

trm12 = mat_times_coeffs(E1t,mtimes_coeffs(C1,fctrD_new));
trm22 = mat_times_coeffs(E1t,mtimes_coeffs(C2,fctrD_new));

trm13 = coeffs_times_mat(mtimes_coeffs(fctrC_new,D1),E1);
trm23 = coeffs_times_mat(mtimes_coeffs(fctrC_new,D2),E1);

P0 = mtimes_coeffs(fctrC0,D0);
P1 = plus_coeffs(plus_coeffs(trm1,trm12),trm13);
P2 = plus_coeffs(plus_coeffs(trm1,trm22),trm23);


P_str = struct();
P_str.C = [P0.C(:)'; P1.C(:)'; P2.C(:)'];
P_str.dim = [m,n];
P_str.dom = dom1;
P_str.deg = 3*d+1;
P_str.vars = [C_str.vars(:,1),D_str.vars(:,2)];

end



%% Extract function-valued 1D PI operators representing the 2D PI operator
function [P0,P1,P2] = get1Dop(P)
% Construct 1D operators P0, P1, P2:L2^{d+1}[a1,b1] --> L2^{d+1}[a1,b1]
% such that
%   (P*x)(s) = (P0(s2)*x(.,s2))(s1)
%               + int_{a2}^{s2} (P1(s2,th2)*(x(.,th2))(s1) dth2
%                   + int_{s2}^{b2} (P2(s2,th2)*x(.,th2))(s1) dth2

dim = P.dim;
deg = P.deg;
vars = P.vars;

P_tmp = struct();
P_tmp.deg = deg;
P_tmp.dim = (deg+1)*dim;
P_tmp.dom = P.dom(2,:);
P_tmp.vars = vars(2,:);

P0 = P_tmp;         P0.C = P.C(1,:);
P1 = P_tmp;         P1.C = P.C(2,:);
P2 = P_tmp;         P2.C = P.C(3,:);

P0.dim = [(deg+1)*dim(1),dim(2)];

end



%% Function for computing sum of operators
function C_str = plus_coeffs(A_str,B_str)
% Compute the sum Cop = Aop+Bop;
C_str = A_str;
for ii=1:numel(C_str.C)
    C_str.C{ii} = A_str.C{ii} + B_str.C{ii};
end
end

%% Function for computing difference of operators
function C_str = minus_coeffs(A_str,B_str)
% Compute the difference Cop = Aop-Bop;
C_str = A_str;
for ii=1:numel(C_str.C)
    C_str.C{ii} = A_str.C{ii} - B_str.C{ii};
end
end


%% Function for multiplying matrix with operator
function C_str = mat_times_coeffs(Amat,B_str)

d = B_str.deg;
m = size(Amat,1);       n = B_str.dim(2);

% Compute the coefficients defining the product
C_str = B_str;
C_str.dim = [m,n];

fctrA = kron(Amat,speye(d+1));
C_str.C{1} = fctrA*B_str.C{1};
C_str.C{2} = fctrA*B_str.C{2};
C_str.C{3} = fctrA*B_str.C{3};

end


function C_str = coeffs_times_mat(A_str,Bmat)

d = A_str.deg;
m = A_str.dim(1);       n = size(Bmat,2);

% Compute the coefficients defining the product
C_str = A_str;
C_str.dim = [m,n];

fctrB = kron(Bmat,speye(d+1));
C_str.C{1} = A_str.C{1}*Bmat;
C_str.C{2} = A_str.C{2}*fctrB;
C_str.C{3} = A_str.C{3}*fctrB;

end
function P_str = mtimes_coeffs(C_str,D_str)
% P_STR = MTIMES_COEFFS(C_STR,D_STR) computes the composition of the
% 3-PI operators defined by structs C_STR and D_STR, such that
%   P_STR = opvar2coeffs(coeffs2opvar(C_STR)*coeffs2opvar(D_STR));
%
% INPUTS
% - C_str,D_str:    1x2 'cell' objects, with each element of both arrays being a
%           d+1 x d+1 array specifying the coefficients defining parameters
%           R.R1 and R.R2 of opvar objects;
%
% OUTPUTS
% - P_str


% Check that coefficients are appropriately specified
dim1 = C_str.dim;   dim2 = D_str.dim;
if dim1(2)~=dim2(1)
    error("Inner dimensions of operators to compose should match.")
end
m = dim1(1);    q = dim1(2);    n = dim2(2);
C0 = C_str.C{1};    C1 = C_str.C{2};      C2 = C_str.C{3};
D0 = D_str.C{1};    D1 = D_str.C{2};      D2 = D_str.C{3};
if any(size(C1)~=size(C2)) || any(size(D1)~=size(D2))
    error("Coefficients representing each parameter should be of same size.")
end
d = size(C1,1)/m-1;
if size(C1,2)~=q*(d+1) || size(D1,1)~=q*(d+1) || size(D1,2)~=n*(d+1)
    error("Coefficients representing the operators should be specified as m*(d+1) x n*(d+1) arrays.")
end

% Check that the domain is appropriately specified
dom1 = C_str.dom;       dom2 = D_str.dom;
if ~isa(dom1,'double') || ~numel(dom1)==2
    error("Spatial domain should be specified as 1x2 array")
end
if any(dom1~=dom2)
    error("Spatial domains of the operators should match.")
end
a = dom1(1);       b = dom1(2);


% Declare the matrices (Iq o A) and (Iq o B)
d_arr = (0:d)+(0:d)'+1;
vals = 1./d_arr(:);
nd = numel(d_arr);
r_idcs = repmat((1:d+1)',[1,d+1]);
r_idcs = r_idcs(:) + repmat((0:q-1)*(d+1),[nd,1]);
c_idcs = repmat((1:d+1),[d+1,1]);
c_idcs = c_idcs(:) + repmat((0:q-1)*(d+1),[nd,1]);
if a==0
    Amat = sparse(q*(d+1),q*(d+1));
elseif a==1
    avals = repmat(vals,[q,1]);
    Amat = sparse(r_idcs,c_idcs,avals,q*(d+1),q*(d+1));
else
    avals = repmat(vals.*(a.^d_arr(:)),[q,1]);
    Amat = sparse(r_idcs,c_idcs,avals,q*(d+1),q*(d+1));
end
if b==0
    Bmat = sparse(q*(d+1),q*(d+1));
elseif b==1
    bvals = repmat(vals,[q,1]);
    Bmat = sparse(r_idcs,c_idcs,bvals,q*(d+1),q*(d+1));
else
    bvals = repmat(vals.*(b.^d_arr(:)),[q,1]);
    Bmat = sparse(r_idcs,c_idcs,bvals,q*(d+1),q*(d+1));
end

% Declare the permutation matrix (Iq o Sd) for Sd such that
%   int Z_{d}(s) Z_{d}(s)' = Sd kron(eye(d+1),Z_{2d+1}(s))
Sd_vals = repmat(vals,[q,1]);
r_idcs = repmat((1:d+1)',[1,d+1]);
r_idcs = r_idcs(:) + repmat((0:q-1)*(d+1),[numel(r_idcs),1]);
c_idcs = (1:d+1)'+(1:d+1)+(2*d+2)*(0:d);
c_idcs = c_idcs(:) + repmat((0:q-1)*(d+1)*(2*d+2),[numel(c_idcs),1]);
Sd = sparse(r_idcs(:),c_idcs(:),Sd_vals(:),q*(d+1),q*(d+1)*(2*d+2));
Sdt = sparse(c_idcs(:),r_idcs(:),Sd_vals(:),q*(d+1)*(2*d+2),q*(d+1));


% Declare the permuation matrix (In o Ed) and (Im o Ed)^T for Ed such that
%   kron(Z_{d}(s),Z_{2d+1}(s)) = Ed*Z_{3d+1}(s)
vals = ones((d+1)*(2*d+2),1);
Ed_vals = repmat(vals(:),[n,1]);
Edt_vals = repmat(vals(:),[m,1]);
r_idcs = (1:(d+1)*(2*d+2))';
c_idcs = (1:2*d+2)'+(0:d);
r_idcs_n = r_idcs(:) + repmat((0:n-1)*(d+1)*(2*d+2),[numel(r_idcs),1]);
c_idcs_n = c_idcs(:) + repmat((0:n-1)*(3*d+2),[numel(c_idcs),1]);
r_idcs_m = c_idcs(:) + repmat((0:m-1)*(3*d+2),[numel(c_idcs),1]);
c_idcs_m = r_idcs(:) + repmat((0:m-1)*(d+1)*(2*d+2),[numel(r_idcs),1]);
Ed = sparse(r_idcs_n(:),c_idcs_n(:),Ed_vals(:),n*(d+1)*(2*d+2),n*(3*d+2));
Edt = sparse(r_idcs_m(:),c_idcs_m(:),Edt_vals(:),m*(3*d+2),m*(d+1)*(2*d+2));

% Declare the permutation matrix (Im o Fd)^T for Fd such that
%   kron(Z_{d}(s),Z_{d}(s)) = Fd*Z_{3d+1}(s)
%Fd_vals = ones(n*(d+1)^2,1);
Fdt_vals = ones(m*(d+1)^2,1);
r_idcs = (1:(d+1)^2);
c_idcs = (1:d+1)' + (0:d);
%r_idcs_n = r_idcs(:) + (0:n-1)*(d+1)^2;
%c_idcs_n = c_idcs(:) + (0:n-1)*(3*d+2);
r_idcs_m = c_idcs(:) + repmat((0:m-1)*(3*d+2),[numel(c_idcs),1]);
c_idcs_m = r_idcs(:) + repmat((0:m-1)*(d+1)^2,[numel(r_idcs),1]);
%Fd = sparse(r_idcs_n(:),c_idcs_n(:),Fd_vals(:),n*(d+1)^2,n*(3*d+2));
Fdt = sparse(r_idcs_m(:),c_idcs_m(:),Fdt_vals(:),m*(3*d+2),m*(d+1)^2);

% Declare the permutation matrix (Iq o Hd) for Hd such that
%   Z_{d}(s)*Zd_{s}' = Hd*(I_{d+1} o Z_{3d+1}(s))
Hd_vals = ones(q*(d+1)^2,1);
r_idcs = repmat((1:d+1)',[1,d+1]);
c_idcs = (1:d+1)' + (0:d)*(3*d+3);
r_idcs_q = r_idcs(:) + (0:q-1)*(d+1);
c_idcs_q = c_idcs(:) + (0:q-1)*(d+1)*(3*d+2);
Hd = sparse(r_idcs_q(:),c_idcs_q(:),Hd_vals(:),q*(d+1),q*(d+1)*(3*d+2));

% Declare the permutation matrix (In o E1) and (Im o E1) such that
%    Z_{d}(s) = E1*Z_{3d+1}(s)
E1_vals = ones(n*(d+1),1);
E1t_vals = ones(m*(d+1),1);
r_idcs_n = (1:d+1)' + (0:n-1)*(d+1);
c_idcs_n = (1:d+1)' + (0:n-1)*(3*d+2);
r_idcs_m = (1:d+1)' + (0:m-1)*(3*d+2);
c_idcs_m = (1:d+1)' + (0:m-1)*(d+1);
E1 = sparse(r_idcs_n(:),c_idcs_n(:),E1_vals,n*(d+1),n*(3*d+2));
E1t = sparse(r_idcs_m(:),c_idcs_m(:),E1t_vals,m*(3*d+2),m*(d+1));
%E1 = kron(eye(n),[speye(d+1),sparse(d+1,2*d+1)]);
%E1t = kron(eye(m),[speye(d+1);sparse(2*d+1,d+1)]);

% Finally, compute the coefficients defining the composition of the
% operators associated with C and D
fctrC0 = Fdt*kron(C0,eye(d+1));
fctrD0 = Hd*kron(D0,eye(3*d+2));
trm1 = E1t*sparse(C2*Bmat*D1-C1*Amat*D2)*E1;
fctrD = Sd*(kron(D2-D1,speye(2*d+2))*Ed);
fctrC = (Edt*kron(C1-C2,speye(2*d+2)))*Sdt;

P0 = fctrC0*D0;
P1 = trm1 + E1t*(C1*(fctrD0+fctrD)) + ((fctrC0+fctrC)*D1)*E1;
P2 = trm1 + E1t*(C2*(fctrD0+fctrD)) + ((fctrC0+fctrC)*D2)*E1;
P = {P0,P1,P2};

P_str = struct();
P_str.C = P;
P_str.dim = [m,n];
P_str.dom = [a,b];
P_str.deg = 3*d+1;
P_str.vars = [C_str.vars(1),D_str.vars(2)];

end


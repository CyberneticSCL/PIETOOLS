function P = mtimes_coeffs(C,D,dom)
% MTIMES_COEFFS Summary of this function goes here
%   Detailed explanation goes here
%
% INPUTS
% - C,D:    1x2 'cell' objects, with each element of both arrays being a
%           d+1 x d+1 array specifying the coefficients defining parameters
%           R.R1 and R.R2 of opvar objects;
% - dom:    1x2 array of type 'double', specifying the spatial interval on
%           which the opvar objects are defined;


% Check that coefficients are appropriately specified
C1 = C{1};      C2 = C{2};
D1 = D{1};      D2 = D{2};
if size(C1,1)~=size(C1,2)
    error("Coefficient matrices should be square.")
elseif any(size(C1)~=size(C2))
    error("Coefficients representing each parameter should be of same size.")
elseif any(size(C1)~=size(D1)) || any(size(C2)~=size(D2))
    error("Coefficients representing the operators should be of the same size.")
end

% Check that the domain is appropriately specified
if ~isa(dom,'double') || ~numel(dom)==2
    error("Spatial domain should be specified as 1x2 array")
end
a = dom(1);       b = dom(2);

d = size(C1,1)-1;

% Declare the matrices A and B
d_arr = (0:d)+(0:d)'+1;
if a==0
    Amat = sparse(d+1,d+1);
elseif a==1
    Amat = (1./d_arr);
else
    Amat = (1./d_arr).*a.^d_arr;
end
if b==0
    Bmat = sparse(d+1,d+1);
elseif b==1
    Bmat = (1./d_arr);
else
    Bmat = (1./d_arr).*b.^d_arr;
end

% Declare the permutation matrix Sd such that
%   int Z_{d}(s) Z_{d}(s)' = Sd kron(eye(d+1),Z_{2d+1}(s))
vals = 1./d_arr;
r_idcs = repmat((1:d+1)',[1,d+1]);
c_idcs = (1:d+1)'+(1:d+1)+(2*d+2)*(0:d);
Sd = sparse(r_idcs(:),c_idcs(:),vals(:),d+1,(d+1)*(2*d+2));
Sdt = sparse(c_idcs(:),r_idcs(:),vals(:),(d+1)*(2*d+2),d+1);


% Declare the permuation matrix Ed such that
%   kron(Z_{d}(s),Z_{2d+1}(s)) = Ed*Z_{3d+1}(s)
vals = ones((d+1)*(2*d+2),1);
r_idcs = (1:(d+1)*(2*d+2))';
c_idcs = (1:2*d+2)'+(0:d);
Ed = sparse(r_idcs(:),c_idcs(:),vals(:),(d+1)*(2*d+2),3*d+2);
Edt = sparse(c_idcs(:),r_idcs(:),vals(:),3*d+2,(d+1)*(2*d+2));

% Declare the permutation matrix E1 such that
%    Z_{d}(s) = E1*Z_{3d+1}(s)
E1 = [speye(d+1),sparse(d+1,2*d+1)];
E1t = [speye(d+1);sparse(2*d+1,d+1)];

% Finally, compute the coefficients defining the composition of the
% operators associated with C and D
trm1 = E1t*sparse(C2*Bmat*D1-C1*Amat*D2)*E1;
fctrD = Sd*(kron(D2-D1,speye(2*d+2))*Ed);
fctrC = (Edt*kron(C1-C2,speye(2*d+2)))*Sdt;

P1 = trm1 + E1t*(C1*fctrD) + (fctrC*D1)*E1;
P2 = trm1 + E1t*(C2*fctrD) + (fctrC*D2)*E1;
P = {P1,P2};


end


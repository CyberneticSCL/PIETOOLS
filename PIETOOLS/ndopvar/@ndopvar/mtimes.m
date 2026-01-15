function Pop = mtimes(Qop,Rop)
% POP = MTIMES(QOP,ROP) returns the 'ndopvar' object representing the
% composition of the PI operators defined by QOP and ROP.
%
% INPUTS
% - Qop:    m x p 'ndopvar' object;
% - Rop:    p x n 'ndopvar' object;
%
% OUTPUTS
% - Pop:    m x n 'ndopvar' object representing the composition of the
%           operators defined by Aop and Bop;
%
% NOTES
% The spatial domains on which the operators defined by Qop and Rop exist
% must match, i.e. Qop.dom=Rop.dom.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes
%
% Copyright (C) 2026 PIETOOLS Team
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
% DJ, 01/15/2026: Initial coding



% Separately deal with the cases of matrix x operator or operator x matrix
% composition
if isa(Qop,'double')
    Pop = mat_times_ndopvar(Qop,Rop);
    return
elseif ~isa(Qop,'ndopvar') && ~isa(Qop,'nopvar')
    error("Only composition of 'ndopvar' objects with 'nopvar' objects is currently supported.")
end
if isa(Rop,'double')
    Pop = ndopvar_times_mat(Qop,Rop);
    return
elseif ~isa(Rop,'ndopvar') && ~isa(Rop,'nopvar')
    error("Only composition of 'ndopvar' objects with 'nopvar' objects is currently supported.")
end
if isa(Qop,'ndopvar') && isa(Rop,'ndopvar')
    error("Composition of two decision variable operators is not supported.")
elseif isa(Rop,'ndopvar')
    % Support for this case will need to be added
    error("Composition of 'nopvar' x 'ndopvar' is currently not supported.")
end


% Check that the inner dimensions of the operators match.
dim1 = Qop.dim;   dim2 = Rop.dim;
if dim1(2)~=dim2(1)
    error("Inner dimensions of operators to compose should match.")
end
m = dim1(1);    p = dim1(2);    n = dim2(2);

if ~isempty(Rop.dvarname)
    error("Composition of known operator with decision variable operator is currently not supported.")
end

% Check that the domain is appropriately specified
dom1 = Qop.dom;       dom2 = Rop.dom;
if ~isa(dom1,'double') || ~numel(dom1)==2
    error("Spatial domain should be specified as 1x2 array")
end
if any(dom1~=dom2)
    error("Spatial domains of the operators should match.")
end
N = size(dom1,1);

% Prohibit cases where the operators are of different monomial degrees
% --> support for such cases will have to be included as well!
deg1 = Qop.deg;     deg2 = Rop.deg;
if any(deg1~=deg2)
    error("Degrees of the opvar2d objects must match.")
end

% Proceed with composition just along the first dimension
d = deg1(1);
a = dom1(1,1);       b = dom1(1,2);


% Declare the matrices (Iq o A) and (Iq o B), where
% B-A = int_{a}^{b} Zd(s)*Zd(s)^T ds 
d_arr = (0:d)+(0:d)'+1;
vals = 1./d_arr(:);
r_idcs = repmat((1:d+1)',[1,d+1]);
c_idcs = repmat((1:d+1),[d+1,1]);
if a==0
    Amat = sparse(p*(d+1),p*(d+1));
elseif a==1
    Amat = spIkron(p,r_idcs,c_idcs,vals,d+1,d+1);
else
    avals = vals.*(a.^d_arr(:));
    Amat = spIkron(p,r_idcs,c_idcs,avals,d+1,d+1);
end
if b==0
    Bmat = sparse(p*(d+1),p*(d+1));
elseif b==1
    Bmat = spIkron(p,r_idcs,c_idcs,vals,d+1,d+1);
else
    bvals = vals.*(b.^d_arr(:));
    Bmat = spIkron(p,r_idcs,c_idcs,bvals,d+1,d+1);
end

% Declare the permutation matrix (Iq o Sd) for Sd such that
%   int Z_{d}(s) Z_{d}(s)' = Sd kron(eye(d+1),Z_{2d+1}(s))
r_idcs = repmat((1:d+1)',[1,d+1]);
c_idcs = (1:d+1)'+(1:d+1)+(2*d+2)*(0:d);
Sd = spIkron(p,r_idcs,c_idcs,vals,d+1,(d+1)*(2*d+2));
Sdt = spIkron(p,c_idcs,r_idcs,vals,(d+1)*(2*d+2),d+1);

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
Hd = spIkron(p,r_idcs(:),c_idcs(:),Hd_vals,(d+1),(d+1)*(3*d+2));

% Declare the permutation matrices (In o E1) and (Im o E1)^T such that
%    Z_{d}(s) = E1*Z_{3d+1}(s)
E1 = spIkron(n,(1:d+1),(1:d+1),1,d+1,3*d+2);
E1t = spIkron(m,(1:d+1),(1:d+1),1,3*d+2,d+1);


% % Deal with the base case: N=1
if N==1
    ndvars = numel(Qop.dvarname);
    if ndvars>0
        % Upscale the matrices to account for the decision variables
        E1t = spIkron(ndvars+1,E1t);
        Fdt = spIkron(ndvars+1,Fdt);
        Edt = spIkron(ndvars+1,Edt);

        % Pre- and post-multiply with appropriate commutation matrices to
        % factor out the decision variables
        Pmat1 = commat(m*(3*d+2),ndvars+1,'transpose');
        Pmat2 = commat(ndvars+1,m*(d+1),'transpose');
        E1t = Pmat1*E1t*Pmat2;
        Fdt = Pmat1*Fdt*kron(Pmat2,speye(d+1));
        Edt = Pmat1*Edt*kron(Pmat2,speye(2*d+2));
    end

    fctrC0 = Fdt*kron(Qop.C{1},eye(d+1));
    fctrD0 = Hd*kron(Rop.C{1},eye(3*d+2));
    trm1 = E1t*sparse(Qop.C{3}*Bmat*Rop.C{2}-Qop.C{2}*Amat*Rop.C{3})*E1;
    fctrD = Sd*(kron(Rop.C{3}-Rop.C{2},speye(2*d+2))*Ed);
    fctrC = (Edt*kron(Qop.C{2}-Qop.C{3},speye(2*d+2)))*Sdt;
    
    P0 = fctrC0*Rop.C{1};
    P1 = trm1 + E1t*(Qop.C{2}*(fctrD0+fctrD)) + ((fctrC0+fctrC)*Rop.C{2})*E1;
    P2 = trm1 + E1t*(Qop.C{3}*(fctrD0+fctrD)) + ((fctrC0+fctrC)*Rop.C{3})*E1;

    Pop = ndopvar();
    Pop.C = {P0;P1;P2};
    Pop.dom = [a,b];
    Pop.deg = 3*d+1;
    Pop.dvarname = Qop.dvarname;
    Pop.vars = [Qop.vars(:,1),Rop.vars(:,2)];
    return
end

% % Otherwise, decompose the ND operator into 3 (N-1)D operators
[C0,C1,C2] = splitNDop(Qop);
[D0,D1,D2] = splitNDop(Rop);

% Perform composition only along the first dimension
C0_I = kronI(C0,d+1);
D0_I = kronI(D0,3*d+2);
C12 = C1-C2;
D21 = D2-D1;
C12_I = kronI(C12,2*d+2);
D21_I = kronI(D21,2*d+2);


% Finally, compute the coefficients defining the composition of the
% operators associated with Qop and Rop
fctrC0 = Fdt*C0_I;
fctrD0 = Hd*D0_I;
CBD = C2*(Bmat*D1);
CAD = C1*(Amat*D2);
AB_diff = CBD-CAD;
trm1 = E1t*(AB_diff*E1);
fctrD = Sd*(D21_I*Ed);
fctrC = (Edt*C12_I)*Sdt;

fctrC_new = fctrC0+fctrC;
fctrD_new = fctrD0+fctrD;

trm12 = E1t*(C1*fctrD_new);
trm22 = E1t*(C2*fctrD_new);

trm13 = (fctrC_new*D1)*E1;
trm23 = (fctrC_new*D2)*E1;

P0 = fctrC0*D0;
P1 = trm1+trm12+trm13;
P2 = trm1+trm22+trm23;


Pop = ndopvar();
sz_C = size(P0.C);
Pop.C = [reshape(P0.C,[1,sz_C]); 
           reshape(P1.C,[1,sz_C]);
           reshape(P2.C,[1,sz_C])];
%P_out.dim = [m,n];
Pop.dom = dom1;
Pop.deg = 3*deg1+1;
Pop.vars = [Qop.vars(:,1),Rop.vars(:,2)];
Pop.dvarname = [Qop.dvarname; Rop.dvarname];

end



%% Extract function-valued (N-1)D PI operators representing the ND operator
function [P0,P1,P2] = splitNDop(P)
% Construct (N-1)D operators,
%   P0, P1, P2:L2^{m*(d1+1)}[a_rem,b_rem] --> L2^{n*(d1+1)}[a_rem,b_rem]
% for [a_rem,b_rem] = [a2,b2] x ... x [aN,bN], such that
%   (P*x)(s) = (P0(s_rem)*x(.,s_rem))(s1)
%               + int_{a1}^{s1} (P1(s1,t1)*(x(t1,.))(s_rem) dt1
%                   + int_{s1}^{b1} (P2(s1,t1)*x(t1,.))(s_rem) dt1,
% for s = (s1,s_rem).

dim = P.dim;
deg = P.deg;
vars = P.vars;

P_tmp = ndopvar();
P_tmp.deg = deg(2:end);
P_tmp.dim = (deg(1)+1)*dim;
P_tmp.dom = P.dom(2:end,:);
P_tmp.vars = vars(2:end,:);
P_tmp.dvarname = P.dvarname;

sz_C = [size(P.C),1];

P0 = P_tmp;         P0.C = reshape(P.C(1,:),sz_C(2:end));
P1 = P_tmp;         P1.C = reshape(P.C(2,:),sz_C(2:end));
P2 = P_tmp;         P2.C = reshape(P.C(3,:),sz_C(2:end));

end





%% Function for multiplying matrix with operator
function Cop = mat_times_ndopvar(Amat,Bop)
% COP = MAT_TIMES_NDOPVAR(AMAT,BOP) returns an 'ndopvar' object
% representing the composition of the matrix AMAT with the operator defined 
% by BOP;

% if ~isempty(Bop.dvarname)
%     error("Composition with decision variable operators is currently not supported.")
% end
q = numel(Bop.dvarname);
d = prod(Bop.deg+1);

% Compute the coefficients defining the product
Cop = Bop;
fctrA = kron(Amat,speye(d*(q+1)));
for ii=1:numel(Cop.C)
    Cop.C{ii} = fctrA*Bop.C{ii};
end

end


%% Function for multiplying operator with matrix
function Cop = ndopvar_times_mat(Aop,Bmat)
% COP = NDOPVAR_TIMES_MAT(AOP,BMAT) returns an 'ndopvar' object
% representing the composition of the operator defined by AOP with the
% matrix BMAT;

if Aop.dim(2)~=size(Bmat,1)
    error("Inner dimensions of operators to compose should match.")
else
    p = size(Bmat,1);
end

% Compute the coefficients defining the product
Cop = Aop;
for ii=1:numel(Aop.C)
    nZ = size(Aop.C{ii},2)/p;
    fctrB = kron(Bmat,speye(nZ));
    Cop.C{ii} = Aop.C{ii}*fctrB;
end

end

function [Tops,Rops,Kops,Fop,Qops] = compute_Tmap_1D(n,E,F,dom,vars,orders)
% [Tops,Rops,Kops,Fop,Qops] = COMPUTE_TMAP_1D(n,E,F,dom,vars,orders)
% takes a matrix E and polynomial F defining BCs on domain dom=[a,b] as
%   sum_{j=1}^{d} E_{i,j}*[D_{j}u_{j}(a); D_{j}u_{j}(b)]
%       + sum_{j=0}^{d} int_{a}^{b} F_{i,j}(s)*u_{j}(s) ds = w_{i},
% for i in {1,...,q}, where D_{j}u = [u; u_{s}; ...; (d/ds)^{j-1}u], and
% u_{j} is a function in W_{2}^{j,n_{j}}[a,b], for state dimensions
% n=(n_{0},...,n_{d}), and returns operators such that, for all
% u=(u_{0},...,u_{d}) satisfying the BCs, with v0=Fop*u and
% v1=(u_{0}; u_{1}'; ...; (d/ds)^{d}u_{d}),
%     u = Tops.T0*v0 + Tops.T1*v1 + Tops.Tw*w,      
%     Kops.K1*v1 = Kops.Kw*w,
% and, for all 0<=i<=j<=d,
%     (d/ds)^{i} u_{j} = Rops.R0{j+1,i+1}*v0 + Rops.R1{j+1,i+1}*v1
%                                                 + Rops.Rw{j+1,i+1}*w,
% and conversely, for all v0 in R^{m} and v1 in L2^{|n|}[a,b] such that
% Kops.K1*v1=Kops.Kw*w, the function u defined by u = Tops.T0*v0 +
% Tops.T1*v1 + Tops.Tw*w satisfies the BCs, with Fop*u=v0, and
% (d/ds)^{j}u_{j}=v1_{j} for all j. Here the "fundamental state" [v0;v1]
% is free of all boundary conditions except Kops.K1*v1=Kops.Kw*w.
%
% INPUT
% - n:  1 x d+1 array specifying the vector dimension of the different
%       state variables, n = (n0,...,nd);
% - E:  q x 2*n_tot array of type 'double', n_tot=sum((1:d).*n(2:end)), for
%       q>=n_tot the number of BCs, defining the contribution of the
%       boundary states to the boundary conditions, with the columns of the
%       block E_{i,j} corresponding to [D_{j}u_{j}(a); D_{j}u_{j}(b)];
% - F:  q x sum(n) array of type 'polynomial' or 'double' specifying
%       the integral term in the BCs. Can depend on at most 1 variable;
% - dom:    1x2 array specifying the lower and upper boundaries [a,b] of
%           the interval on which u(s) is defined;
% - vars:   (optional) 1x2 array of type 'polynomial', specifying the
%           spatial variable s and a corresponding dummy variable theta;
% - orders: (optional) array of the orders i of the derivatives (d/ds)^{i}
%           for which to compute Rops and Qops.Rp, 0:d (all) by default;
%
% OUTPUT
% All operators are 'copvar' objects, each a 1 x 1 container holding a
% single 'sopvar' block, so that every operator maps between just one
% input and one output space (R^{p} or L2^{p}[a,b]). The spatial variable
% is vars(1) (or the variable of F, or the default of 'opvar'), and the
% dummy variable is <vars(1)>_dum, as 'sopvar' does not record the dummy.
% Here m = n_tot-rank(G) is the number of auxiliary values, mK = q-rank(G)
% the number of constraints on the fundamental state (mK=m if q=n_tot),
% and |n| = sum(n).
% - Tops:   struct with fields
%               T0: R^{m} --> L2^{|n|}, the map v0 --> u (empty if m=0),
%               T1: L2^{|n|} --> L2^{|n|}, the map v1 --> u,
%               Tw: R^{q} --> L2^{|n|}, the map w --> u;
% - Rops:   struct with fields R0, R1, Rw, each a (d+1) x (d+1) cell, with
%               R0{j+1,i+1}: R^{m} --> L2^{n_{j}}  (empty if m=0),
%               R1{j+1,i+1}: L2^{|n|} --> L2^{n_{j}},
%               Rw{j+1,i+1}: R^{q} --> L2^{n_{j}},
%           mapping v0, v1, and w to (d/ds)^{i} u_{j}, for 0<=i<=j<=d.
%           Left empty if n_{j}=0 or i>j;
% - Kops:   struct with fields
%               K1: L2^{|n|} --> R^{mK}, v1 --> int_{a}^{b} K(s)*v1(s) ds,
%               Kw: R^{q} --> R^{mK}, w --> J2*w,
%           such that K1*v1=Kw*w for the fundamental state of u. The
%           constraint does not depend on v0. Both are empty if mK=0;
% - Fop:    L2^{|n|} --> R^{m}, u --> v0 = int_{a}^{b} F3(s)*u(s) ds. Empty
%           if m=0. Not a struct, as v0 does not depend on w. The choice of
%           F3 is not unique. If m=|n|, F3 is the mean of the state,
%           F3 = I/(b-a), provided this defines a valid map, and otherwise
%           F3 = P'*Phi, as described below;
% - Qops:   (optional) struct splitting the map w --> u through the
%           Taylor coefficients c=Cw*w at s=a of the polynomial part of u
%           due to w:
%               w --> c --> u,      Tops.Tw = Tp*Cw,
%                                   Rops.Rw{j+1,i+1} = Rp{j+1,i+1}*Cw.
%           The entries of c are the boundary values Lambda_a u =
%           [D_{1}u_{1}(a);...;D_{d}u_{d}(a)] of single derivatives of
%           single components, instead of the arbitrary combinations in
%           the rows of w. Fields:
%               Cw: n_tot x q array of type 'double', with c = Cw*w;
%               Tp: R^{n_tot} --> L2^{|n|}, c --> ebar(s-a)*c, with ebar
%                   the Taylor basis of the differentiable components;
%               Rp: (d+1) x (d+1) cell, with Rp{j+1,i+1}: R^{n_tot} -->
%                   L2^{n_{j}} the map c --> (d/ds)^{i} u_{j} for the
%                   polynomial ebar(s-a)*c;
%           Left empty if n_{j}=0 or i>j, as for Rops;
%
% NOTES
% The BCs are first converted to canonical form (compute_canonical_BCs),
%   G*[D_{1}u_{1}(a);...;D_{d}u_{d}(a)] + int_{a}^{b} H(s)*v1(s) ds = w,
% where the first r=rank(G) rows (G1,H1) are such that G1 has full row
% rank, and the last mK rows constrain only v1, defining K. The remaining
% boundary values are fixed by m auxiliary values v0 = int_{a}^{b} F3(s)*
% u(s) ds, with F3 chosen such that Ghat = [G1; int F3(s)*ebar(s-a) ds] is
% invertible: F3 = I/(b-a) if m=|n| and possible, and F3 = P'*Phi
% otherwise, where P selects the m columns of G1 that are not pivots. Then
%   T0 = ebar*C0,       T1 = ebar*C1,               Tw = ebar*Cw,
%   C0 = Ghat\[0;I],    C1(t) = -Ghat\[H1;H3](t),   Cw = Ghat\[J1;0],
% with ebar the Taylor basis of the differentiable components. The 'sopvar'
% blocks are constructed directly from the polynomial parameters, as
%   (P*x)(s) = (I kron Z_L(s)) * C0 * (I kron Z_R) x                (R0)
%              + int_{a}^{s} (I kron Z_L(s)) * C1 * (I kron Z_R(t)) x(t) dt
%              + int_{s}^{b} (I kron Z_L(s)) * C2 * (I kron Z_R(t)) x(t) dt
% for monomial vectors Z_L, Z_R, and sparse coefficient matrices C0,C1,C2.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - compute_Tmap_1D
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
% DJ, 11/02/2026: Initial coding;
% DJ, 11/06/2026: Allow q>=n_tot BCs. The number m=n_tot-r of auxiliary
%   values (v0) and the number mK=q-r of constraints on v1 differ for
%   q>n_tot, so K1,Kw are declared whenever mK>0 (was m>0), and w has q
%   entries (was n_tot);


% % Convert BCs to canonical form (which checks the inputs),
% %   G*Lambda_a u + int_{a}^{b} H(s)*v1(s) ds = w
% % with J*[BCs] reordered such that the first r rows involve Lambda_a u.
if nargin<5
    vars = [];
end
[G,H,p,J,vars] = compute_canonical_BCs(n,E,F,dom,vars);
% A 'sopvar' does not record the dummy variable, so the convention
% <var1>_dum is used (as in opvar2sopvar), whatever dummy was specified
var1 = vars(1);     var2 = pvar([char(var1.varname{1}),'_dum']);
a = dom(1);         b = dom(2);

d = numel(n)-1;
if nargin<6
    orders = 0:d;
end
n_tot = sum((1:d).*n(2:end));       % number of boundary values
n_all = sum(n);                     % dimension of the state [u_0;...;u_d]
r = numel(p);
m = n_tot-r;
% The q-r rows of G beyond its rank r are constraints on v1, outnumbering
% the m=n_tot-r auxiliary values (the boundary values not fixed by the BCs)
% if q>n_tot
q = size(E,1);                                                              % DJ, 11/06/2026
mK = q-r;
c = setdiff(1:n_tot,p);             % non-pivot columns of G

G1 = G(1:r,:);          H1 = H(1:r,:);          J1 = J(1:r,:);
K = H(r+1:end,:);       J2 = J(r+1:end,:);

% Selection matrices for the components of the state, 
%   u_j = S{j+1}*[u_0;...;u_d],
% and for their boundary values, 
%   Lambda_a u_j = R{j}*Lambda_a u, j>=1
idx_G0 = cumsum([0,(1:d).*n(2:end)]);
idx_F0 = cumsum([0,n]);
S = cell(1,d+1);        R = cell(1,d);
for j=0:d
    S{j+1} = zeros(n(j+1),n_all);
    S{j+1}(:,idx_F0(j+1)+(1:n(j+1))) = eye(n(j+1));
    if j>=1
        R{j} = zeros(j*n(j+1),n_tot);
        R{j}(:,idx_G0(j)+(1:j*n(j+1))) = eye(j*n(j+1));
    end
end


% % % Construct the coordinates, and the auxiliary values
% % Declare the polynomial ebar(s-a), mapping the boundary values to the
% % state,
% %     ebar(z) = [0; diag(e_{1,1}(z),...,e_{d,1}(z))],
% % where e_{j,1}(z) = [p_0(z)*I, ..., p_{j-1}(z)*I] for p_k(z) = z^k/k!
ebar = polynomial(zeros(n_all,n_tot));
for j=1:d
    for k=0:j-1
        ek = zeros(1,j);        ek(k+1) = 1;
        ebar = ebar + (var1-a)^k/factorial(k)*(S{j+1}'*kron(ek,eye(n(j+1)))*R{j});
    end
end

% % Declare the polynomial Phi(s) such that
% % int_{a}^{b} Phi(s)*ebar(s-a) ds = I,
% %     Phi(s) = [0, diag(Phi_1(s),...,Phi_d(s))],
% % with Phi_j(s) = col(phi_0(s)*I,...,phi_{j-1}(s)*I). Here,
% % Here phi_k is dual to p_l(s-a), l=0,...,j-1, defined using the inverse
% % of the the Hilbert matrix.
Phi = polynomial(zeros(n_tot,n_all));
for j=1:d
    x = (var1-a)/(b-a);
    for k=0:j-1
        phi_k = 0;
        for l=0:j-1
            c_kl = (-1)^(k+l)*(k+l+1)*nchoosek(j+k,j-1-l)*nchoosek(j+l,j-1-k)*nchoosek(k+l,k)^2;
            phi_k = phi_k + c_kl*x^l;
        end
        phi_k = factorial(k)/(b-a)^(k+1)*phi_k;
        ek = zeros(j,1);        ek(k+1) = 1;
        Phi = Phi + phi_k*(R{j}'*kron(ek,eye(n(j+1)))*S{j+1});
    end
end

% % Auxiliary values v0 = int_{a}^{b} F3(s)*u(s) ds. Any F3 for which
% %     Ghat = [G1;G3],     G3 = int_{a}^{b} F3(s)*ebar(s-a) ds,
% % is invertible will do. By default, we take
% %     F3 = P'*Phi,
% % where P=I(:,c) selects the non-pivot columns. Then G3 = P' 
% % (as int_{a}^{b} Phi(s)*ebar(s-a) ds = I), so that, as G1(:,p) = I, 
% % for Pi=I(:,[p,c]),
% %     Ghat*Pi = [I, G1(:,c); 0, I],
% % which can be inverted in closed form.
Pi = eye(n_tot);
Pi = Pi(:,[p,c]);
Ghat_inv = Pi*[eye(r),-G1(:,c);zeros(m,r),eye(m)];
if m>0
    % If m=|n|, first try the mean of the state, F3 = I/(b-a), unless Ghat
    % is not invertible for it (e.g. if u_0 has no boundary values to
    % recover), and revert to the default
    F3_opts = {Phi(c,:)};
    if m==n_all
        F3_opts = [{eye(n_all)/(b-a)},F3_opts];
    end
    for k=1:numel(F3_opts)
        % G3 and H3 such that v0 = G3*Lambda_a u + int H3(s)*v1(s) ds
        F3 = F3_opts{k};
        [G3,H3] = compute_canonical_BCs(n,zeros(m,2*n_tot),F3,dom,vars,'nosplit');
        Ghat_mean = [G1;G3];
        if k<numel(F3_opts) && rcond(Ghat_mean)>1e-8
            Ghat_inv = inv(Ghat_mean);
            break
        end
    end
end


% % % Compute the maps v0, v1, w --> u, where u = T0*v0 + T1*v1 + Tw*w

% Coefficients of the boundary values,
%   Lambda_a u = C0*v0 + Cw*w + int C1(s)*v1(s) ds
C0 = Ghat_inv(:,r+1:end);
Cw = Ghat_inv(:,1:r)*J1;
C1 = polynomial(zeros(n_tot,n_all));
if r>0
    C1 = C1 -Ghat_inv(:,1:r)*H1;
end
if m>0
    C1 = C1 -Ghat_inv(:,r+1:end)*H3;
end
C1 = subs(C1,var1,var2);

% u = ebar(s-a)*Lambda_a u + P0*v1(s)
%       + int_{a}^{s} Pbar(s-theta)*v1(theta) dtheta
T0 = ebar*C0;           Tw = ebar*Cw;       T1 = ebar*C1;
P0 = blkdiag(eye(n(1)),zeros(n_all-n(1)));
Pbar = polynomial(zeros(n_all,n_all));
for j=1:d
    Pbar = Pbar + (var1-var2)^(j-1)/factorial(j-1)*(S{j+1}'*S{j+1});
end


% % % Declare the operators as 'copvar' objects, directly from their
% % % polynomial parameters. Each is a map between a single pair of spaces,
% % % R^{p} or L2^{p}[a,b], so a single 'sopvar' block suffices.
vname = char(var1.varname{1});
mk = @(polys,q,pp,in_L2,out_L2) copvar({make_sopvar(polys,q,pp,vname,[a,b],in_L2,out_L2)});

% u = T0*v0 + T1*v1 + Tw*w. T1 has multiplier P0 and kernels
% Pbar(s-theta)+T1(s,theta) on [a,s], T1(s,theta) on [s,b]
Tops = struct('T0',[],'T1',[],'Tw',[]);
Tops.T1 = mk({P0,Pbar+T1,T1},n_all,n_all,true,true);
Tops.Tw = mk({Tw},n_all,q,false,true);
if m>0
    Tops.T0 = mk({T0},n_all,m,false,true);
end

% K1*v1 = Kw*w, with K1*v1 = int K(s)*v1(s) ds, and Kw = J2
Kops = struct('K1',[],'Kw',[]);
Fop = [];
% The constraints exist whenever mK>0, also if there are no auxiliary
% values (m=0)
if mK>0                                                                     % DJ, 11/06/2026
    Kops.K1 = mk({K},mK,n_all,true,false);
    Kops.Kw = mk({J2},mK,q,false,false);
end
if m>0
    Fop = mk({F3},m,n_all,true,false);
end

% (d/ds)^{i} u_j = (d/ds)^{i} S_j*[T0*v0 + T1*v1 + Tw*w], where for i<j,
%   (d/ds)^{i} int_{a}^{s} p_{j-1}(s-theta)*v1_j(theta) dtheta
%       = int_{a}^{s} p_{j-1-i}(s-theta)*v1_j(theta) dtheta,
% and for i=j, (d/ds)^{j} u_j = S_j*v1.
Rops = struct('R0',{cell(d+1,d+1)},'R1',{cell(d+1,d+1)},'Rw',{cell(d+1,d+1)});
for j=0:d
    nj = n(j+1);
    if nj==0
        continue
    end
    for i=orders(orders<=j)
        if i==j
            R0ji = zeros(nj,m);     Rwji = zeros(nj,q);
            R1ji = {S{j+1},zeros(nj,n_all),zeros(nj,n_all)};
        else
            T0ji = S{j+1}*T0;       
            T1ji = S{j+1}*T1;       
            Twji = S{j+1}*Tw;
            for k=1:i
                T0ji = diff(T0ji,var1);     
                T1ji = diff(T1ji,var1);     
                Twji = diff(Twji,var1);
            end
            R0ji = T0ji;            
            Rwji = Twji;
            R1ji = {0*S{j+1},(var1-var2)^(j-1-i)/factorial(j-1-i)*S{j+1}+T1ji,T1ji};
        end
        if m>0
            Rops.R0{j+1,i+1} = mk({R0ji},nj,m,false,true);
        end
        Rops.R1{j+1,i+1} = mk(R1ji,nj,n_all,true,true);
        Rops.Rw{j+1,i+1} = mk({Rwji},nj,q,false,true);
    end
end


% % % Split Tw=ebar*Cw and Rw into the Taylor coefficients c=Cw*w and the
% % % maps c --> (d/ds)^{i} u_{j} for the polynomial ebar(s-a)*c, so that
% % % the coefficients can be treated per component
if nargout>=5
    Qops = struct('Cw',Cw,'Tp',mk({ebar},n_all,n_tot,false,true),'Rp',{cell(d+1,d+1)});
    for j=0:d
        nj = n(j+1);
        if nj==0
            continue
        end
        for i=orders(orders<=j)
            Rpji = zeros(nj,n_tot);         % (d/ds)^{j} of a polynomial of degree <j
            if i<j
                Rpji = S{j+1}*ebar;
                for k=1:i
                    Rpji = diff(Rpji,var1);
                end
            end
            Qops.Rp{j+1,i+1} = mk({Rpji},nj,n_tot,false,true);
        end
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function P = make_sopvar(polys,q,p,vname,dom,in_L2,out_L2)
% P = MAKE_SOPVAR(POLYS,Q,P,VNAME,DOM,IN_L2,OUT_L2) declares the 'sopvar'
% mapping L2^{p}[dom] (if IN_L2) or R^{p} to L2^{q}[dom] (if OUT_L2) or
% R^{q}.
%
% INPUT
% - polys:  cell of q x p parameters, of type 'polynomial' or 'double':
%           - L2-->L2:  {R0(s),R1(s,theta),R2(s,theta)}, defining
%                R0(s)*x(s) + int_{a}^{s} R1(s,theta)*x(theta) dtheta
%                           + int_{s}^{b} R2(s,theta)*x(theta) dtheta,
%                       in variables s=vname and theta=<vname>_dum;
%           - R-->L2:   {Q2(s)}, defining the map x --> Q2(s)*x;
%           - L2-->R:   {Q1(s)}, defining the map
%                       x --> int_{a}^{b} Q1(s)*x(s) ds
%                       (the polynomial is in the variable named vname);
%           - R-->R:    {M}, a q x p matrix of type 'double';
% - q,p:    integers, the output and input dimensions of the operator;
% - vname:  char, the name of the spatial variable;
% - dom:    1 x 2 array, the interval [a,b] of the variable;
% - in_L2,out_L2:  logical, true if the input (output) space is
%           L2^{p}[dom] (L2^{q}[dom]), and false if it is R^{p} (R^{q});
%
% OUTPUT
% - P:      'sopvar' object, with the operator specified by polys;
%
% NOTES
% The parameters of 'sopvar' are sparse matrices C such that, for monomial
% exponents ZL and ZR, a kernel K(s,theta) = sum_{k,l} C(i*nZL+k,j*nZR+l)*
% s^ZL(k)*theta^ZR(l) in entry (i,j) (zero-based i,j; components outer,
% monomials inner). A multiplier R0(s) is stored with theta-degree 0 only,
% as required by the canonical multiplier form of 'sopvar'.
dname = [vname,'_dum'];
vars = struct;      vars.in = {};           vars.out = {};
dom_s = struct;     dom_s.in = zeros(0,2);  dom_s.out = zeros(0,2);
ZL = cell(1,0);     ZR = cell(1,0);
if ~in_L2 && ~out_L2
    params = {full(double(polys{1}))};
else
    if in_L2
        vars.in = {vname};      dom_s.in = dom;
    end
    if out_L2
        vars.out = {vname};     dom_s.out = dom;
    end
    % Extract the terms of all parameters (q x p polynomials in vname, and
    % for L2-->L2 also in dname): for the nt terms of parameter k, the
    % exponents dgL{k} (nt x 1) in the output variable and dgR{k} in the
    % input variable, and the matrix coefficients cf{k} (nt x q*p), with
    % row t the column-wise vectorized coefficient of term t
    nP = numel(polys);
    dgL = cell(1,nP);   dgR = cell(1,nP);   cf = cell(1,nP);
    for k=1:nP
        Pk = polynomial(polys{k});
        if ~isequal(size(Pk),[q,p])
            error("Parameter has dimension %dx%d, but %dx%d was expected.",size(Pk,1),size(Pk,2),q,p)
        end
        cf{k} = full(Pk.coefficient);       dm = full(Pk.degmat);
        dgL{k} = zeros(size(cf{k},1),1);    dgR{k} = dgL{k};
        for v=1:numel(Pk.varname)
            if strcmp(Pk.varname{v},vname)
                if out_L2
                    dgL{k} = dm(:,v);
                else
                    dgR{k} = dm(:,v);
                end
            elseif strcmp(Pk.varname{v},dname) && in_L2 && out_L2
                dgR{k} = dm(:,v);
            else
                error("Parameter depends on unexpected variable '"+Pk.varname{v}+"'.")
            end
        end
    end
    % Common monomial bases. A multiplier has theta-degree 0, so 0 is
    % included in ZR for 3-PI operators.
    dL = unique([0;vertcat(dgL{:})]);
    dR = unique([0;vertcat(dgR{:})]);
    if ~out_L2
        dL = 0;
    end
    if ~in_L2
        dR = 0;
    end
    nZL = numel(dL);    nZR = numel(dR);
    if out_L2
        ZL = {dL};
    end
    if in_L2
        ZR = {dR};
    end
    % Assemble the coefficient matrices
    params = cell(1,nP);
    for k=1:nP
        [~,li] = ismember(dgL{k},dL);
        [~,ri] = ismember(dgR{k},dR);
        [tt,ij,vals] = find(cf{k});
        ii = mod(ij-1,q)+1;     jj = floor((ij-1)/q)+1;
        params{k} = sparse((ii-1)*nZL+li(tt),(jj-1)*nZR+ri(tt),vals,q*nZL,p*nZR);
    end
end
P = sopvar(params,vars,ZL,ZR,dom_s,[q,p]);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

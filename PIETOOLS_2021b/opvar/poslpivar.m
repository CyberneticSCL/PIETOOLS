function [prog,Pop] = poslpivar(prog,n,I,d,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop] = poslpivar(prog,n,I,d,options) declares 
% a positive 4-PI decision operator Pop, with fields:
%
% P = Q11 int(gs,s,a,b)
% Q(s) = g(s) Q12 Z1(s) + int( gth Q13 Z2(th,s) dth, s, b) + int( gth Q14 Z3(th,s) dth, a, s)
% R0(s) = gs Z1(s) Q22 Z1(s)
% R1(s,th) = gs Z1(s) Q23 Z2(s,th) + gth Z3(th,s) Q42 Z1(th) + 
%            int(geta Z2(eta,s)' Q33 Z2(eta,th) deta,s,b)+int(geta Z3(eta,s) Q43 Z2(eta,th) deta,th,s)
%           +int(geta Z3(eta,s) Q44 Z3(eta,th) deta,a,th)
% R2(s,th) = R1(th,s)',
%
% where
%
% Q = [ Q_{11}  Q_{12} Q_{13} Q_{14}]
%     [ Q_{21}  Q_{22} Q_{23} Q_{24}] >0
%     [ Q_{31}  Q_{32} Q_{33} Q_{34}]
%     [ Q_{41}  Q_{42} Q_{43} Q_{44}]
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
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - poslpivar
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding MMP, SS, DJ  - 09/26/2021

% % % Set-up % % %

% % Extract the input arguments
switch nargin
    case 2
        error('Not enough inputs!')
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

% Specify the degrees of the monomials
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

% Specify the spatial domain (s,th \in [I(1),I(2)]
if length(I)~=2
    error('I must be a length 2 vector')
end

% Extract the size of the object: P\in R^(n1xn1), R0\in L_2^(n2xn2)
n = reshape(n,[1,2]);
n1 = n(1);  n2 = n(2);
if n1==0
    if n2==0
        error('Error in posopvar: All dimensions are zero')
    end
end


% % To reduce complexity, allow certain terms to be excluded

% excludeL is a length-4 binary vector of terms to exclude
excludeL = options.exclude;

% In separable case R1=R2, so Qi3 = Qi4 and Q33 = Q44, so we may exclude
% Qi4 entirely, and use Qi3 instead.
if options.sep ==1 
    excludeL(4)=1;
end

% If there is no finite-dimensional component or infinite-dimensional
% components, we may exclude the associated term
if n1==0
    excludeL(1)=1;
end
if n2==0
    excludeL(2:4)=[1 1 1];
end

if all(excludeL)
    error('You''re creating an empty dopvar! Please change options.exclude and/or options.sep and try again.')
end

% % Define the variable of multiplier function

% Define the primary variables (R0 is function of s, R1 and R2 of s,theta)
pvar s theta;
var1 = s; var2 = theta;

% Define a dummy variable for integration (replaces eta)
sss =  polynomial(1,1,{'sss'},[1 1]);

% Define the multiplier function to be used later
if options.psatz==1 % use a function that is positive only on domain I
    gs = (var1-I(1))*(I(2)-var1);
    gth = (var2-I(1))*(I(2)-var2);
    geta = (sss-I(1))*(I(2)-sss);
else % use a function that is positive everywhere
    gs = polynomial(1);
    gth = polynomial(1);
    geta = polynomial(1);
end

% % % Construct the monomial bases Z and positive matrix Q % % %

% % Build the monomials
% Construct Z1(s)
nZ1 = d{1}+1;
Z1degmat = (0:d{1})';   % maximal degree d{1}
Z1coeff = speye(nZ1);
Z1varname = var1.varname;
Z1matdim = [nZ1 1];
Z1s = polynomial(Z1coeff,Z1degmat,Z1varname,Z1matdim);

% Construct Z2(s,th) and Z3(s,th)
% In this implementation, Z2 will have degree d{2,2} in th and degree 
% d{2,1} in s, and max degree of s+th is d{2,3}. Similarly for Z3(s,th)
Z2degmat = [repmat((0:d{2}(1))',d{2}(2)+1,1),vec(repmat(0:d{2}(2),d{2}(1)+1,1))];
Z2degmat(sum(Z2degmat,2)>d{2}(3),:)= [];
nZ2=size(Z2degmat,1);
Z2coeff = speye(nZ2);
Z2varname = [var1.varname; var2.varname];
Z2matdim = [nZ2 1];
Z2sth=polynomial(Z2coeff,Z2degmat,Z2varname,Z2matdim);

Z3degmat = [repmat((0:d{3}(1))',d{3}(2)+1,1),vec(repmat(0:d{3}(2),d{3}(1)+1,1))];
Z3degmat(sum(Z3degmat,2)>d{3}(3),:)= [];
nZ3=size(Z3degmat,1);
Z3coeff = speye(nZ3);
Z3varname = [var1.varname; var2.varname];
Z3matdim = [nZ3 1];
Z3sth=polynomial(Z3coeff,Z3degmat,Z3varname,Z3matdim);

% % Build the associated block matrices
% bZ1 = [Z1,0 ,...,0 ;]   etc.
%       [0 ,Z1,...,0 ;]
%       [: ,: , . ,: ;]
%       [0 ,0 ,...,Z1 ]

% The size of the matrices must be such that bZ1'*Q22*bZ1 \in L_2^(n2xn2) 
%nBZ1 = n2*nZ1;
%nBZ2 = n2*nZ2;
%nBZ3 = n2*nZ3;

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

% Swap the variables s and th for later purposes
%bZ2ths = var_swap(bZ2sth,var1,var2);
%bZ3ths = var_swap(bZ3sth,var1,var2);


% % Declare the positive matrix variable

% We are going to compute bZL'*Q*bZR
includeL = ~excludeL;
ZL = cell(1,sum(includeL));
ZR = cell(1,sum(includeL));

Z1th = subs(Z1s,var1,var2);
Z2etath = subs(Z2sth,var1,sss);
Z2etas = subs(Z2etath,var2,var1);
Z3etath = subs(Z3sth,var1,sss);
Z3etas = subs(Z3etath,var2,var1);

mdim = [];
ndim = [];
indx = 1;
if includeL(1)
    ZL{indx} = 1;
    ZR{indx} = 1;
    mdim = [mdim;n1];
    ndim = [ndim;n1];
    indx = indx+1;
end
if includeL(2)
    ZL{indx} = Z1s;
    ZR{indx} = Z1th;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(3)
    ZL{indx} = Z2etas;
    ZR{indx} = Z2etath;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(4)
    ZL{indx} = Z3etas;
    ZR{indx} = Z3etath;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
end
[prog,N] = sosquadvar(prog,ZL,ZR,mdim,ndim,'pos');

include_indx = cumsum(includeL);
    


% % This option is currently not available
if options.diag
    disp('Poslpivar option ''diag'' is currently not available')
end
% This option is really only intended for time-delay systems, wherein we
% know these terms are zero and can be omitted
% if options.diag==1 && options.exclude(2)~=1 % if the diagonal override is used, make $Q22$ block diagonal
%     Q22temp=Q{2,2}; %this is of dimension nBZ1=n2*nZth=K*diag*nZth
%     K=n2/diag; % this is the number of blocks
%     for i=1:K
%         irange=(diag*(i-1)*nZs+1):(i*diag*nZs);    % n2/
%         Q22temp(irange,irange)=zeros(diag*nZs); % parts of Q22 which are not Zero
%     end
%     for i=1:K
%         irange=(diag*(i-1)*nZs+1):(i*diag*nZs);    % n2/
%         for j=1:K
%             if i~=j
%                 jrange=(diag*(j-1)*nZs+1):(j*diag*nZs);    % n2/
%                 Q{2,2}(irange,jrange)=zeros(length(irange),length(jrange));
%                 % eliminate these terms so fewer variables appear in the polynomial manipulations
%             end
%         end
%     end
%     prog=sosmateq(prog,Q22temp); %constrain the off-diagonal terms to be zero
% end


% % % Build the positive dopvar object % % %

% Since the object must be self-adjoint, we need only define Pop.P = P,
% Pop.Q1 = QT, Pop.R.R0 = R0 and Pop.R.R1 = R1
P = dpvar(zeros(n1));
QT = dpvar(zeros(n1,n2));
R0 = dpvar(zeros(n2));
R1 = R0;

% Define P, Q1 = QT, and Q2 = QT'
if excludeL(1)==0
    %P = P+Q{1,1}*int(gs,var1,I(1),I(2));
    %P = P + int(gs * ZQZconstruct(poly1,poly1,Q{1,1},[n1,n1],var1,var2,sss,1,0),var1,I(1),I(2));
    i1 = include_indx(1);     j1 = include_indx(1);
    P = P + N{i1,j1} * int(gs,var1,I(1),I(2));
    if excludeL(2)==0
        %QT = QT + Q{1,2}*gs*bZ1s;
        %QT = QT + gs * ZQZconstruct(poly1,Z1s,Q{1,2},[n1,n2],var1,var2,sss,0,0);
        i1 = include_indx(1);     j1 = include_indx(2);
        QT = QT + gs * subs(N{i1,j1},var2,var1);
    end
    if excludeL(3)==0 && ~options.sep==1
        %QT = QT + int(Q{1,3}*gth*bZ2ths,var2,var1,I(2));
        %QT = QT + int(gth * ZQZconstruct(poly1,Z2sth,Q{1,3},[n1,n2],var1,var2,sss,0,0),var2,var1,I(2));
        i1 = include_indx(1);     j1 = include_indx(3);
        QT = QT + int(geta * subs(N{i1,j1},var2,var1),sss,var1,I(2));
    elseif excludeL(3)==0 && options.sep==1
        %QT = QT + int(Q{1,3}*gth*bZ2ths,var2,I(1),I(2));
        %QT = QT + int(gth * ZQZconstruct(poly1,Z2sth,Q{1,3},[n1,n2],var1,var2,sss,0,0),var2,I(1),I(2));
        i1 = include_indx(1);     j1 = include_indx(3);
        QT = QT + int(geta * subs(N{i1,j1},var2,var1),sss,I(1),I(2));
    end
    if excludeL(4)==0
        %QT = QT + int(Q{1,4}*gth*bZ3ths,var2,I(1),var1);
        %QT = QT + int(gth * ZQZconstruct(poly1,Z2sth,Q{1,4},[n1,n2],var1,var2,sss,0,0),var2,I(1),var1);
        i1 = include_indx(1);     j1 = include_indx(4);
        QT = QT + int(geta * subs(N{i1,j1},var2,var1),sss,I(1),var1);
    end
end

% Define terms in R0 and R1 related to the first monomial basis
if excludeL(2)==0
    %R0 = gs * bZ1s'*Q{2,2}*bZ1s;
    %R0 = gs * ZQZconstruct(Z1s,Z1s,Q{2,2},[n2,n2],var1,var2,sss,1,0);
    i1 = include_indx(2);     j1 = include_indx(2);
    R0 = R0 + gs * subs(N{i1,j1},var2,var1);
    if excludeL(3)==0 && ~options.sep==1
        %R1 = R1 + gs * bZ1s'*Q{2,3}*bZ2sth;
        %R1 = R1 + gs * ZQZconstruct(Z1s,Z2sth,Q{2,3},[n2,n2],var1,var2,sss,0,0);
        i1 = include_indx(2);     j1 = include_indx(3);
        R1 = R1 + gs * subs(N{i1,j1},sss,var1);
    elseif excludeL(3)==0 && options.sep==1
        %R1 = R1 + gs * bZ1s'*Q{2,3}*bZ2sth + gth * bZ2ths'*Q{2,3}.'*bZ1th;
        %R1 = R1 + gs * ZQZconstruct(Z1s,Z2sth,Q{2,3},[n2,n2],var1,var2,sss,0,0) ...
        %        + gth * ZQZconstruct(Z1s,Z2sth,Q{2,3},[n2,n2],var1,var2,sss,0,1);
        i1 = include_indx(2);     j1 = include_indx(3);
        R1 = R1 + gs * subs(N{i1,j1},sss,var1) + gth * subs(N{j1,i1},sss,var2);
    end
    if excludeL(4)==0
        %R1 = R1 + gth * bZ3ths'*Q{2,4}.'*bZ1th;
        %R1 = R1 + gth * ZQZconstruct(Z1s,Z2sth,Q{2,4},[n2,n2],var1,var2,sss,0,1);
        i1 = include_indx(4);     j1 = include_indx(2);
        R1 = R1 + gth * subs(N{i1,j1},sss,var2);
    end
end

% Define terms in R1 related to the second monomial basis
if excludeL(3)==0 && ~options.sep==1
        %R1 = R1 + int(geta * bZ2etas'*Q{3,3}*bZ2etath,sss,var1,I(2)); 
        %R1 = R1 + int(geta * ZQZconstruct(Z2sth,Z2sth,Q{3,3},[n2,n2],var1,var2,sss,1,0),sss,var1,I(2));
        i1 = include_indx(3);     j1 = include_indx(3);
        R1 = R1 + int(geta * N{i1,j1},sss,var1,I(2));
    if excludeL(4)==0
        %R1 = R1 + int(geta * bZ3etas'*Q{3,4}.'*bZ2etath,sss,var2,var1); 
        %R1 = R1 + int(geta * ZQZconstruct(Z2sth,Z3sth,Q{3,4},[n2,n2],var1,var2,sss,0,1),sss,var2,var1);
        i1 = include_indx(4);     j1 = include_indx(3);
        R1 = R1 + int(geta * N{i1,j1},sss,var2,var1);
    end
elseif excludeL(3)==0 && options.sep==1
    %R1 = R1 + int(geta * bZ2etas'*Q{3,3}*bZ2etath,sss,I(1),I(2));
    %R1 = R1 + int(geta * ZQZconstruct(Z2sth,Z2sth,Q{3,3},[n2,n2],var1,var2,sss,1,0),sss,I(1),I(2));
    i1 = include_indx(3);     j1 = include_indx(3);
    R1 = R1 + int(geta * N{i1,j1},sss,I(1),I(2));
end

% Define terms in R1 related to the third monomial basis
if excludeL(4)==0
    %R1 = R1 + int(geta * bZ3etas'*Q{4,4}*bZ3etath,sss,I(1),var2);
    %R1 = R1 + int(geta * ZQZconstruct(Z3sth,Z3sth,Q{4,4},[n2,n2],var1,var2,sss,1,0),sss,I(1),var2);
    i1 = include_indx(4);     j1 = include_indx(4);
    R1 = R1 + int(geta * N{i1,j1},sss,I(1),var2);
end

% Define R2
if options.sep==1
    R2 = R1;
else
    R2 = varswap(R1,var1,var2).';
end

% Construct the positive dpvar Pop
dopvar Pop;
Pop.P = P;
Pop.Q1 = QT;
Pop.Q2 = QT.'; Pop.R.R0 = R0; Pop.R.R1 = R1; Pop.R.R2 = R2;

Pop.dim = [n',n'];
Pop.I = I;
Pop.var1 = var1;
Pop.var2 = var2;

end
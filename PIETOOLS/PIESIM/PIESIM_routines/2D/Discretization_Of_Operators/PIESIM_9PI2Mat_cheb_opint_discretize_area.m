%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_9PI2Mat_cheb_opint_discretize_area.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of
% double-integrated blocks (R11, R12, R21, R22) of 9PI operator in 2D
%
% Called by 'PIESIM_9PI2Mat_cheb_2D.m' 
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial block (R11, R12, R21 or R22) of 9PI operator in 2D
% corresponding to a single solution state
% rsize - the number of rows in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the solution state on the left-hand side of the ODE matrix system corresponding to the discrete block in question)
% csize - the number of columns in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the right-hand side of the ODE matrix system corresponding to the discrete block in question)
% lim  - limits of 2D integraton (for example, [-1 sx;-1 sy]) 
%
% Outputs:
% A - Chebyshev
% discretizaiton of R11, R12, R21 or R22 block that represents a block of a
% square total matrix operator for time-advancement of the
% spatially-discretized PIE solution (square ODE system matrix)
% A_nonsquare - Chebyshev
% discretizaiton of R11, R12, R21 or R22 block that represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution (nonsquare transformation matrix)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024

function [Afull, Afull_nonsquare]=PIESIM_9PI2Mat_cheb_opint_discretize_area(N, R, rsize, csize, lim)
pvar s1 s1_dum s2 s2_dum;

if isa(R,'polynomial')
    deg=R.maxdeg;
else
    deg=1;
end

Reval=zeros(deg+2,deg+2,deg+2,deg+2);

chebgrid=cos(pi*(0:deg+1)/(deg+1));

if isa(R,'polynomial')
if ismember('s2_dum',R.varname)
    Ry=subs(R,s2_dum,chebgrid);
else
    Ry=R*ones(1,deg+2);
end
if ismember('s2',R.varname)
    Rynu=subs(Ry,s2,chebgrid);
else
    Rynu=repmat(Ry,deg+2,1);
end
if ismember('s1_dum',R.varname) & ismember('s1',R.varname)
    for j=1:deg+2
        for i=1:deg+2
        Rmat=subs(subs(Rynu(i,j),s1_dum,chebgrid),s1,chebgrid);
        Reval(:,:,i,j)=Rmat(:,:);
        end
    end
elseif ismember('s1',R.varname)
   for j=1:deg+2
        for i=1:deg+2
        Re=subs(Rynu(i,j),s1,chebgrid);
        Reval(:,:,i,j)=repmat(Re',1,deg+2);
        end
   end
   else
       for j=1:deg+2
        for i=1:deg+2
       Re=subs(Rynu(i,j),s1_dum,chebgrid);
       Reval(:,:,i,j)=repmat(Re,deg+2,1);
        end
       end
   end
else
    Reval=R*ones(deg+2,deg+2,deg+2,deg+2);
end

Reval=double(Reval);
Rcheb=fcgltran4d(Reval);

for q=1:deg+2
    for p=1:deg+2
        A(1:N+1,1:csize,p,q)=PIESIM_integral_projection1D(N, Rcheb(:,:,p,q), csize, lim(1,:));
    end
end

 for k=1:csize
     for i=1:N+1
         Rycheb(:,:)=A(i,k,:,:);
         A2D(i,k,1:N+1,1:csize)=PIESIM_integral_projection1D(N, Rycheb, csize, lim(2,:));
     end
 end



acheb_all=permute(A2D,[1 3 2 4]);
acheb=acheb_all(1:rsize,1:rsize,:,:);
Afull_nonsquare=reshape(acheb_all,(N+1)^2,csize^2);
Afull=reshape(acheb,rsize^2,csize^2);





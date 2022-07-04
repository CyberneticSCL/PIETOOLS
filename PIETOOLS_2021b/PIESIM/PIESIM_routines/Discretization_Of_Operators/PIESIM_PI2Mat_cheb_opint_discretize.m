%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opint_discretize.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation for the integrative Q1 operator  
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% nx - number of ODE states
% Rop - Q1 operator of dimension nx x ns
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for PDE, see Peet & Peet 2021 paper
%
% Outputs:
% A - discretization matrix of the Q1 operator: 
% dimension nx x n0(N+1) n1N n2(N-1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_30_2022
function A=PIESIM_PI2Mat_cheb_opint_discretize(N, nx, Rop, p);
pvar s

ns=size(p,2);

if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

Reval=zeros(nx,ns,deg+2);
acheb=zeros(nx,ns,deg+2);

chebgrid=cos(pi*(0:deg+1)/(deg+1));
if isa(Rop,'polynomial')
for j=1:ns
    for i=1:nx
    Reval(i,j,:)=(subs(Rop(i,j),s,chebgrid))';
    acheb(i,j,:)=fcht(double(Reval(i,j,:)));
    end 
end 
else
for j=1:ns
for i=1:nx
    Reval(i,j,:)=Rop(i,j);
    acheb(i,j,:)=fcht(double(Reval(i,j,:)));
end
end
end


  
  for j=1:ns
        Norder=N-p(j);
       
  A=zeros(nx,Norder);
  
        for k=0:Norder
            for i=1:nx
 vecint=0;
 for m=0:deg+1
% Integrate Rop (T_(k+m) and T_|k-m|) polynomials on [-1,1]

index=[k+m,abs(k-m)];
for kk=index
 if (kk==0)
     vecint=vecint+0.5*acheb(i,j,m+1)*(1)^1-0.5*acheb(i,j,m+1)*(-1)^1;
 elseif (kk==1)
     vecint=vecint+0.5*0.25*acheb(i,j,m+1)*((1)^0+(1)^2)-0.5*0.25*acheb(i,j,m+1)*((-1)^0+(-1)^2);
 else
    vecint=vecint+0.5*0.5*acheb(i,j,m+1)/(kk+1)*(1)^(kk+1)-0.5*0.5*acheb(i,j,m+1)/(kk-1)*(1)^(kk-1)...
    -0.5*0.5*acheb(i,j,m+1)/(kk+1)*(-1)^(kk+1)+0.5*0.5*acheb(i,j,m+1)/(kk-1)*(-1)^(kk-1);
end % if
end % for kk=index
end % for m
       
        A(i,k+1)=vecint;
            end
        end
        
        A_cell{j}=double(A);
  end
        
        A=cell2mat(A_cell);

  




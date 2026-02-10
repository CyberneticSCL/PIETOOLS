<<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_PI2Mat_cheb_opint_discretize.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opint_discretize.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation for the integrative Q1 operator  
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop - Q1 operator of dimension no x ns
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for PDE, see Peet & Peet 2021 paper
%
% Outputs:
% A - discretization matrix of the Q1 operator: 
% dimension no x n0(N+1) n1N n2(N-1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_30_2022% 
% YP -  added support for an arbitary variable name - 04_16_2024

function A=PIESIM_PI2Mat_cheb_opint_discretize(N, Rop, p);
pvar var s

ns=length(p);
no=size(Rop,1);

if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

Reval=zeros(no,ns,deg+2);
acheb=zeros(no,ns,deg+2);

chebgrid=cos(pi*(0:deg+1)/(deg+1));
if isa(Rop,'polynomial')
for j=1:ns
    for i=1:no
        if (isempty(Rop.varname))
            var=s;
        else
            var=Rop.varname;
        end
    Reval(i,j,:)=(subs(Rop(i,j),var,chebgrid))';
    acheb(i,j,:)=fcht(double(Reval(i,j,:)));
    end 
end 
else
for j=1:ns
for i=1:no
    Reval(i,j,:)=Rop(i,j);
    acheb(i,j,:)=fcht(double(Reval(i,j,:)));
end
end
end


  
  for j=1:ns
        Norder=N-p(j);
       
  A=zeros(no,Norder);
  
        for k=0:Norder
            for i=1:no
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

  



========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_opint_cheb.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representation for the integrative Q1 operator  
% Discretizes a FULL OPERATOR, but does not provide a PIE-2-PDE component
% (since it maps to ODE states)
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop - Q1 operator of dimension no x ns
% p - vector of dimension 1xns -
% a "degree of smoothness" structure for PDE, see Peet & Peet 2021 paper
%
% Outputs:
% A - discretization matrix of the Q1 operator: 
% dimension no x n0(N+1) n1N n2(N-1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_30_2022% 
% YP -  added support for an arbitary variable name - 04_16_2024

function A=PIESIM_PI2Mat_opint_cheb(N, Rop, p);
pvar var s

ns=length(p);
no=size(Rop,1);

if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

Reval=zeros(no,ns,deg+2);
acheb=zeros(no,ns,deg+2);

chebgrid=cos(pi*(0:deg+1)/(deg+1));
if isa(Rop,'polynomial')
for j=1:ns
    for i=1:no
        if (isempty(Rop.varname))
            var=s;
        else
            var=Rop.varname;
        end
    Reval(i,j,:)=(subs(Rop(i,j),var,chebgrid))';
    acheb(i,j,:)=fcht(double(Reval(i,j,:)));
    end 
end 
else
for j=1:ns
for i=1:no
    Reval(i,j,:)=Rop(i,j);
    acheb(i,j,:)=fcht(double(Reval(i,j,:)));
end
end
end


  
  for j=1:ns
        Norder=N-p(j);
       
  A=zeros(no,Norder);
  
        for k=0:Norder
            for i=1:no
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

  



>>>>>>>> Stashed changes:PIETOOLS/PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_PI2Mat_opint_cheb.m

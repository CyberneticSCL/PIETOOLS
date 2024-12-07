%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI2Mat_cheb_opint_discretize.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - a dicrete matrix representation of an integrative portion
% of the 3PI operator
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R1 - integrative operator R1 of size 1x1
% R2 - integrative operator R2 of size 1x1
% p - scalar - a "degree of smoothness" for the function on which
% R1, R2 operators act
% var1 - primary spatial variable used in the integral operator defined by
%           R1 and R2;
% var2 - dummy variable (for integration) used in the integral operator
%           defined by R1 and R2;
%
% Output:
% A - a non-square matrix of size (N+1)x(N-p+1) that represents a Chebyshev
% discretizaiton of the integrative portion of the PI operator
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_28_2022
% DJ, 12/07/2024: Update to avoid hardcoded polynomial variable;
function A=PIESIM_3PI2Mat_cheb_opint_discretize(N, R1, R2, p, var1, var2)

s = var1;                                                                   % DJ, 12/07/2024
theta = var2;

Norder=N-p;

A(1:N+1,1:Norder+1)=0;

if isa(R1,'polynomial')
    deg1=R1.maxdeg;
else
    deg1=1;
end

if isa(R2,'polynomial')
    deg2=R2.maxdeg;
else
    deg2=1;
end

maxdeg=max(deg1,deg2);

chebgrid=cos(pi*(0:maxdeg+1)/(maxdeg+1));
if isa(R1,'polynomial')
    switch(R1.nvars)
        case 2
        % Two variables present
            Reval1=subs(subs(R1,theta,chebgrid),s,chebgrid);
        case 1
        % One variable present
            var=cell2mat(R1.varname);
            if (var==s.varname{1})
                 Re1=subs(R1,s,chebgrid);
                 Reval1=repmat(Re1',1,maxdeg+2);
            else
                Re1=subs(R1,theta,chebgrid);
                Reval1=repmat(Re1,maxdeg+2,1);
            end
        case 0 
            % Scalar polynomial
            Reval1=R1*ones(maxdeg+2);
    end
else
    Reval1=R1*ones(maxdeg+2);
end

Reval1=double(Reval1);
Rcheb1=fcgltran2d(Reval1,1);

if isa(R2,'polynomial')
switch(R2.nvars)
        case 2
        % Two variables present
            Reval2=subs(subs(R2,theta,chebgrid),s,chebgrid);
        case 1
        % One variable present
            var=cell2mat(R2.varname);
            if (var==s.varname{1})
                 Re2=subs(R2,s,chebgrid);
                 Reval2=repmat(Re2',1,maxdeg+2);
            else
                Re2=subs(R2,theta,chebgrid);
                Reval2=repmat(Re2,maxdeg+2,1);
            end
         case 0 
         % Scalar polynomial
         Reval2=R2*ones(maxdeg+2);
end
else
    Reval2=R2*ones(maxdeg+2);
end
Reval2=double(Reval2);
Rcheb2=fcgltran2d(Reval2,1);


for k=0:Norder
    % Finding a column vector of Chebyshev coefficient for the column k+1
    % of matrix A
vecmult1=zeros(N+1,1);
vecmult2=zeros(N+1,1);
for j=0:maxdeg+1
    for m=0:maxdeg+1
% Integrate R1 (T_(k+m) and T_|k-m|) polynomials on [-1,s]
% Integrate R2 (T_(k+m) and T_|k-m|) polynomials on [s,1]

vecint1=zeros(N+1,1);
vecint2=zeros(N+1,1);

i=[k+m,abs(k-m)];
for kk=i
 if (kk==0)
     vecint1(2)=vecint1(2)+0.5*Rcheb1(j+1,m+1);
     vecint2(2)=vecint2(2)-0.5*Rcheb2(j+1,m+1);
     vecint1(1)=vecint1(1)-0.5*Rcheb1(j+1,m+1)*(-1)^1;
     vecint2(1)=vecint2(1)+0.5*Rcheb2(j+1,m+1)*(1)^1;

  elseif (kk==1)
     vecint1(1)=vecint1(1)+0.5*0.25*Rcheb1(j+1,m+1);
     vecint1(3)=vecint1(3)+0.5*0.25*Rcheb1(j+1,m+1);
     vecint2(1)=vecint2(1)-0.5*0.25*Rcheb2(j+1,m+1);
     vecint2(3)=vecint2(3)-0.5*0.25*Rcheb2(j+1,m+1);
     vecint1(1)=vecint1(1)-0.5*0.25*Rcheb1(j+1,m+1)*((-1)^0+(-1)^2);
     vecint2(1)=vecint2(1)+0.5*0.25*Rcheb2(j+1,m+1)*((1)^0+(1)^2);
else
    if (kk+1+1<=N+1)
    vecint1(kk+1+1)=vecint1(kk+1+1)+0.5*0.5*Rcheb1(j+1,m+1)/(kk+1);
    vecint2(kk+1+1)=vecint2(kk+1+1)-0.5*0.5*Rcheb2(j+1,m+1)/(kk+1);
    end
    if (kk-1+1<=N+1)
    vecint1(kk-1+1)=vecint1(kk-1+1)-0.5*0.5*Rcheb1(j+1,m+1)/(kk-1);
    vecint2(kk-1+1)=vecint2(kk-1+1)+0.5*0.5*Rcheb2(j+1,m+1)/(kk-1);
    end
    vecint1(1)=vecint1(1)-0.5*0.5*Rcheb1(j+1,m+1)/(kk+1)*(-1)^(kk+1)+0.5*0.5*Rcheb1(j+1,m+1)/(kk-1)*(-1)^(kk-1);
    vecint2(1)=vecint2(1)+0.5*0.5*Rcheb2(j+1,m+1)/(kk+1)*(1)^(kk+1)-0.5*0.5*Rcheb2(j+1,m+1)/(kk-1)*(1)^(kk-1);
end % if
end % for kk=i

% Multiply by T_j(s)
   if (j==0) 
       vecmult1=vecmult1+vecint1;
       vecmult2=vecmult2+vecint2;
   else 
       for i=1:N+1
           iorder=i-1;
       if (iorder+j+1<=N+1)
       vecmult1(iorder+j+1)=vecmult1(iorder+j+1)+0.5*vecint1(i);
       vecmult2(iorder+j+1)=vecmult2(iorder+j+1)+0.5*vecint2(i);
       end
       if (abs(iorder-j)+1<=N+1)
       vecmult1(abs(iorder-j)+1)=vecmult1(abs(iorder-j)+1)+0.5*vecint1(i);
       vecmult2(abs(iorder-j)+1)=vecmult2(abs(iorder-j)+1)+0.5*vecint2(i);
       end
       end % i
    end % j
    end % m
end % j
A(:,k+1)=vecmult1(:)+vecmult2(:);
end % k

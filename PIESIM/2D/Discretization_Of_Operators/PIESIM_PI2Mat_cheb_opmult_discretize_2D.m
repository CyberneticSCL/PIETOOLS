%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI2Mat_cheb_opmult_discretize_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A, A_nonsquare - discrete matrix representations of a polymonial
% multiplicative block (R00) of 9PI operator in 2D
%
% Called by 'PIESIM_9PI2Mat_cheb_2D.m' and 'PIESIM_9PI2Mat_cheb_opint_discretize_line.m'
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% R -  polymonial multiplicative block (R00) of 9PI operator in 2D
% corresponding to a single solution state
% rsize - the number of rows in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the solution state on the left-hand side of the ODE matrix system corresponding to the discrete block in question)
% csize - the number of columns in the resulting A matrix block (number of the PIE state Chebyshev coefficients per dimension in the right-hand side of the ODE matrix system corresponding to the discrete block in question)
%
% Outputs:
% A - Chebyshev
% discretizaiton of a R00 block that represents a block of a
% square total matrix operator for time-advancement of the
% spatially-discretized PIE solution (square ODE system matrix)
% A_nonsquare - Chebyshev
% discretizaiton of a R00 block that represents a block of a
% nonsquare total matrix operator for reconstruction of the primary (PDE)
% solution (nonsquare transformation matrix)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 4_16_2024

function [A, A_nonsquare]=PIESIM_PI2Mat_cheb_opmult_discretize_2D(N, R, rsize, csize)

pvar s1 s2;

chebgrid=cos(pi*(0:csize-1)/(csize-1));


    if isa(R,'polynomial')
    Reval=subs(subs(R,s1,chebgrid)',s2,chebgrid);
    else
    Reval=R*ones(1,csize^2);
    end

    acheb=reshape(fcgltran2d(double(Reval),1),[],1);

A(1:rsize^2,1:csize^2)=0;
A_nonsquare(1:(N+1)^2,1:csize^2)=0;

for i=1:csize^2
    id_x=mod(i-1,csize);
    id_y=floor((i-id_x+1)/csize);
    for j=1:length(acheb)
        jd_x=mod(j-1,csize);
        jd_y=floor((j-jd_x+1)/csize);
            pos=(jd_y+id_y)*csize+(jd_x+id_x)+1;
            if (pos<=(N+1)^2)
            A_nonsquare(pos,i)=A_nonsquare(pos,i)+0.25*acheb(j);
            end
            pos=(abs(jd_y-id_y))*csize+(jd_x+id_x)+1;
            if (pos<=(N+1)^2)
            A_nonsquare(pos,i)=A_nonsquare(pos,i)+0.25*acheb(j);
            end
            pos=(jd_y+id_y)*csize+abs(jd_x-id_x)+1;
            if (pos<=(N+1)^2)
            A_nonsquare(pos,i)=A_nonsquare(pos,i)+0.25*acheb(j);
            end
            pos=(abs(jd_y-id_y))*csize+abs(jd_x-id_x)+1;
            if (pos<=(N+1)^2)
            A_nonsquare(pos,i)=A_nonsquare(pos,i)+0.25*acheb(j);
            end
        end 
        end 

for i=1:csize^2
    id_x=mod(i-1,csize);
    id_y=floor((i-id_x+1)/csize);
    for j=1:length(acheb)
        jd_x=mod(j-1,csize);
        jd_y=floor((j-jd_x+1)/csize);
            pos=(jd_y+id_y)*csize+(jd_x+id_x)+1;
            if (pos<=rsize^2) 
            A(pos,i)=A(pos,i)+0.25*acheb(j);
            end
            pos=(abs(jd_y-id_y))*csize+(jd_x+id_x)+1;
            if (pos<=rsize^2) 
            A(pos,i)=A(pos,i)+0.25*acheb(j);
            end
            pos=(jd_y+id_y)*csize+abs(jd_x-id_x)+1;
            if (pos<=rsize^2) 
            A(pos,i)=A(pos,i)+0.25*acheb(j);
            end
            pos=(abs(jd_y-id_y))*csize+abs(jd_x-id_x)+1;
            if (pos<=rsize^2) 
            A(pos,i)=A(pos,i)+0.25*acheb(j);
            end
    end
end


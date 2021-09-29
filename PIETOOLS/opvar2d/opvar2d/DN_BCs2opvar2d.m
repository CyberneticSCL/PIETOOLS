function [Ebb] = DN_BCs2opvar2d(Bindx,n11,n22,dom)
% % Construct opvar2d object Ebb describing Dirichlet and Neumann bc's, 
% 0 = Ebb*ubf;

% Conditions are enforce on hypercube described by dom=[a(1),b(1);a(2),b(2)];

% Bindx describes conditions on first order state of size n1
% Bindx must be binary column vector of length 4
% Bindx(i)==1 enforces Dirichlet condition at boundary i \in {1,2,3,4}
% Bindx(1)+Bindx(3) must be equal to 1, enforcing one condition in dimension 2
% Bindx(2)+Bindx(4) must be equal to 1, enforcing one condition in dimension 1


% Dindx describes conditions on second order state of size n2
% Dindx must be binary matrix of size [4,2]
% Dindx(i,1)==1 enforces Dirichlet condition at boundary i \in {1,2,3,4}
% Dindx(i,2)==1 enforces Neumann condition at boundary i \in {1,2,3,4}
% Dindx(1,1)+Dindx(1,2)+Dindx(3,1)+Dindx(3,2) must be equal to 2, enforcing two conditions in dimension 2
% Dindx(2,1)+Dindx(2,2)+Dindx(4,1)+Dindx(4,2) must be equal to 2, enforcing two conditions in dimension 2

% Edges (i) and corners are numbered as follows:
%   
%               b(2)
%
%       2         3         3
%         o---------------o
%         |               |
%         |               |
% a(1)  2 |               | 4  b(1)
%         |               |
%         |               |
%         o---------------o
%       1         1         4
%               
%               a(2)

% This numbering deviates from order in full boundary vector u_bf
% We use the following vector to translate order of corners in code to
% order in u_bf:
bf_order = [1;3;4;2];

Bindx11 = Bindx(1:4*n11,1);
Bindx22 = Bindx(4*n11+1:end,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n11
    %if any(Bindx11(:,2))
    %    error('Neumann type BCs are not allowed for once differentiable states')
    %end
    if sum(Bindx11([1,3]))~=1
        error('Insufficient or too many boundary conditions provided on second dimension')
    end
    if sum(Bindx11([2,4]))~=1
        error('Insufficient or too many boundary conditions provided on first dimension')
    end
    
    B11 = opvar2d();
    B11.I = dom;
    B11.dim = [n11,4*n11;
        n11,2*n11;
        n11,2*n11;
        0,0];
    
    dif = [Bindx11(2:end);Bindx11(1)] - Bindx11(:);
    b0 = (dif==0 & Bindx11==1);
    b0 = bf_order(b0);
    b1 = Bindx11([1,3])==1;
    b2 = Bindx11([2,4])==1;
    
    B11.R00(1:n11,(b0-1)*n11+1:b0*n11) = eye(n11);
    B11.Rxx{1}(1:n11,(b1-1)*n11+1:b1*n11) = eye(n11);
    B11.Ryy{1}(1:n11,(b2-1)*n11+1:b2*n11) = eye(n11);
    
    B11.dim = B11.dim;
else
    B11 = opvar2d();
    B11.I = dom;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Second order conditions
if n22
    if sum(sum(Bindx22([1,3],:)))~=2
        error('Insufficient or too many boundary conditions provided on second dimension')
    end
    if sum(sum(Bindx22([2,4],:)))~=2
        error('Insufficient or too many boundary conditions provided on first dimension')
    end
    
    B22 = opvar2d();
    B22.I = dom;
    B22.dim = [4*n22,16*n22;
        2*n22,4*n22;
        2*n22,4*n22;
        0,0];
    
    if all(Bindx22(:,1))
        dindx = [1;
            8 + 2;
            8 + 3;
            4 + 4];
    elseif all(Bindx22(:,2))
        error('All Neumann conditions... solution not uniquely defined')
    else
        dif1 = [Bindx22(2:end,1);Bindx22(1,1)] - Bindx22(:,1);
        dif2 = [Bindx22(2:end,2);Bindx22(1,2)] - Bindx22(:,2);
        dif22 = [dif2(2:end);dif2(1)] - dif2;
        d00 = dif1==1; % Corners at which to impose u(.,.)=0
        d10 = find(Bindx22([1,3],1)==1 | (Bindx22([1,3],2)==1 & dif22([1,3])==0)); % Corners at which to impose d_1 u(.,.)=0
        d10 = 1 + 2*(d10-1);
        d01 = find(Bindx22([2,4],1)==1 | (Bindx22([2,4],2)==1 & dif22([2,4])==0));  % Corners at which to impose d_2 u(.,.)=0
        d01 = 2 + 2*(d01-1);
        d11 = (dif22==-1 & dif2==0 & Bindx22(:,2)==1); % Corners at which to impose d_1 d_2 u(.,.)=0
        
        bf_order = [1;3;4;2];   % Order of corners in bf vector
        d00 = bf_order(d00);
        d10 = 4 + bf_order(d10);
        d01 = 8 + bf_order(d01);
        d11 = 12 + bf_order(d11);
        
        dindx = [d00;d10;d01;d11];
    end
    
    B22.R00(1:n22,(dindx(1)-1)*n22+1:dindx(1)*n22) = eye(n22);
    B22.R00(n22+1:2*n22,(dindx(2)-1)*n22+1:dindx(2)*n22) = eye(n22);
    B22.R00(2*n22+1:3*n22,(dindx(3)-1)*n22+1:dindx(3)*n22) = eye(n22);
    B22.R00(3*n22+1:4*n22,(dindx(4)-1)*n22+1:dindx(4)*n22) = eye(n22);
    
    dindx_x = [Bindx22(1,1);Bindx22(3,1);Bindx22(1,2);Bindx22(3,2)];
    dindx_x = find(dindx_x==1);
    dindx_y = [Bindx22(2,1);Bindx22(4,1);Bindx22(2,2);Bindx22(4,2)];
    dindx_y = find(dindx_y==1);
    
    B22.Rxx{1}(1:n22,(dindx_x(1)-1)*n22+1:dindx_x(1)*n22) = eye(n22);
    B22.Rxx{1}(n22+1:2*n22,(dindx_x(2)-1)*n22+1:dindx_x(2)*n22) = eye(n22);
    
    B22.Ryy{1}(1:n22,(dindx_y(1)-1)*n22+1:dindx_y(1)*n22) = eye(n22);
    B22.Ryy{1}(n22+1:2*n22,(dindx_y(2)-1)*n22+1:dindx_y(2)*n22) = eye(n22);
    
    
    B22.dim = B22.dim;
    
else
    B22 = opvar2d();
    B22.I = dom;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Combine the objects
Ebb = blkdiag(B11,B22);

end
function [prog,Q]=sosposmatrvar(prog,n,d,vars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SOSPOSMATRVAR(prog,n,d,vars)
%
% This program declares a symbolic positive scalar semidefinite matrix P of 
% size nxn which is positive semidefinite for all values of the variables.
% The matrix has the form
% [Z   ]^T   [Z   ]
% [ Z  ]   Q [ Z  ]
% [  Z ]     [  Z ]
% [   Z]     [   Z], Q>0
%
%Where $Z$ is the vector of monomials in vaiables vars of degree d/2 or
%less. Q is a positive semidefinite matrix.
%
% INPUTS: 
%
% prog - The SOS program to which to attach the variable
% n - dimension of matrix variable
% d - degree of polynomial entries of the matrix
% vars - variables in the entries of the polynomial
%
%
%
% NOTES:
% Distributed with DelayTOOLS
% Compatable with MULTIPOLY and SOSTOOLS as of June 2013
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

% version .03   M. Peet, matthew.peet@inria.fr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MatrixTools - sosposmatrvar
% Release Notes:
% v.03 - compatability with sosposmatr
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ~even(d)
%   error(['Degree of the SOS Polynomial Matrix should be even']);
% end
if ~isvalid(vars)
  error(['vars must be a polynomial variable']);
end

Z=monomials(vars,0:ceil(d/2));
nZ=length(Z);
[prog,P]=sosposmatr(prog,nZ*n);

% I need to creat a map from the variable name index to the position in the matrix
% P. So that the i,j th element of the matrix LLL is the variable with
% index ind(i,j). Since I know how P is constructed, this shouldn't be
% too bad. Assuming I NEVER change how sosposmatr_p works...



% In this routine, we manually calculate M=bZ1th.'*P*bZ1th;
% All variables appear 
avarname=[P.varname; Z.varname];
% In the final polynomial, there will still be e.g. 7260 terms. 

% For the next part, we just need the degmat part of the block-diagonal
% monomial matrix

bZdegmat=repmat(Z.degmat,n,1);

% fortunately, the index of the elements of P correspond to the power of x
% by which they multiply
[Ilist, Jlist]=find(triu(ones(nZ*n)));
% so we just add those monomials to each term
adegmat = [P.degmat bZdegmat(Ilist,:)+bZdegmat(Jlist,:)];


% but now we have to modify the coefficient matrix. We can group this by
% blocks. For a 2x2 matrix, all the coefficients from block 1 will be [1 0 0 0], all the
% coefficients from block 2 will be [0 1 0 0], all the coefficients from
% block 3 will be [0 0 0 1],

[PIlist, PJlist]=find(P.coefficient); % these will all be ones, but PIlist and Pjlist will be relatively short. 
% Positions in Pjlist from 1:nZ get reassigned to 1, nZ+1:2nZ -> 2, etc
row=mod(PJlist-1,nZ*n)+1; % row number
col=ceil(PJlist/(nZ*n)); % column number
newrow=ceil(row/nZ);
newcol=ceil(col/nZ);
newidx=(newcol-1)*n+newrow;
coeff=sparse(PIlist,newidx,1);
amatdim=[n n];
Q=polynomial(coeff,adegmat,avarname,amatdim);

% 
% for jj=1:n %col count
%     for ii=1:n % row count
%         % This is the sub-indexing matrix for the ii,jj block
%         indiijj = indt((ii-1)*nZ+1:(ii)*nZ,(jj-1)*nZ+1:(jj)*nZ);
%         [Isublist, Jsublist]=find(indiijj);
%         Isubslist=Isublist+(ii-1)*nZ;
%         Jsubslist=Jsublist+(jj-1)*nZ;%these are the positions of each term
%         % We need to assign a coefficient to each monomial in degmat. The
%         % problem is that there are lots of terms in degmat.
%     end
% end
% 
% 
% 
% amatdim=[n n];
% 
% 
% for jj=1:n %col count
%     for ii=1:n % row count
%         % form A(ii,jj). This is a scalar formed from Z^T *P((ii-1)*nZ+1:(ii)nZ,(jj-1)*nZ+1:jj*nZ)*Z
%         % to do this indirectly, we will use the sub-indexing matrix
%         indiijj = ind((ii-1)*nZ1+1:(ii)*nZ1,(jj-1)*nZ1+1:(jj)*nZ1);
%         % for each ii,jj, there will be (nZ*(nZ+1))/2 monomials
%         acoeff=2*ones(,1);
%         for jjj=1:nZ1 % column count
%             for iii=1:jjj % row count
%                 count=count+1;
%                 % The ZZZth degmat 
%                 adegmat(count,nvars)=iii+jjj; % this is the degree of Zth(iii)*Zth(jjj)
%                 % if iii=jjj, the coefficient is 1 instead of two
%                 adegmat(count,indiijj(iii,jjj))=1;
%                 if iii=jjj
%                     acoeff(count)=1;
%                 end
%                 % if iii=jjj, the coefficient is 2
%             end
%         end
%                 % there are two coefficients for each term
%                 adegmat(count,indiijj(jjj,iii))=1; % when iii=jjj, we just do it twice
%         A(ii,jj)=polynomial(acoeff,adegmat,avarname,amatdim);
%         % form B(ii,jj)
%         % form CA(ii,jj)
%     end
% end

function prob = Sedumi2Mosek(A,b,c,K)
% Sedumi2Mosek --- Convert Sedumi inputs to the MOSEK prob format.
%
% prob = Sedumi2Mosek(A,b,c,K)
%
% Description: The Sedumi SDP is defined as
% min c^T x :
% Ax=b
% x_i \in K_i
% The elements of K are: K.f=free, K.l=PO, K.q=Lor, K.r=R Lor, K.s=SDP
% We will assume only K.f, K.l and K.s are non-empty. quadratic constraints
% will return an error.
%
% The mosek SDP format is
% min c^Tx + \sum_i <C_i,X_i>
% lbi \le a^Tx + \sum_j <A_(ij),X_j> \le ubi,
%
% which is all contained in the prob structure.
%
% Inputs:
% 
% A: An m x n matrix
% b: A length m column vector
% c: A length n column vector
% K: is the cone structure from the Sedumi input format.
%
% prob: is the input structure accepted by MOSEK
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

% 04/24/14 - MP - initial coding

% A Sedumi to Mosek converter

if isfield(K,'q')
    error('K.q and K.r are supposed to be empty')
end

if ~isfield(K,'f')
    K.f=[];
end
if ~isfield(K,'l')
    K.l=[];
end
nf=sum(K.f); nl=sum(K.l);

Kindex = nf+nl+cumsum([0 (K.s).^2])+1; % vector with the starting indices of the SDP variables
n=size(A,2);
m=size(A,1);

nsdpvar=length(K.s);


lpA=A(:,1:(Kindex(1)-1));
lpc=c(1:(Kindex(1)-1));

% The structure is stored in prob
prob.barc.subj=[];
prob.barc.subk=[];
prob.barc.subl=[];
prob.barc.val=[];

prob.bara.subi=[];
prob.bara.subj=[];
prob.bara.subk=[];
prob.bara.subl=[];
prob.bara.val=[];
for j=1:nsdpvar
    temp=mat(c(Kindex(j):(Kindex(j+1)-1)),K.s(j));
    [I,J,V] = find(tril(temp+temp')/2);
    prob.barc.subj=[prob.barc.subj j*ones(1,length(I))];
    prob.barc.subk=[prob.barc.subk I' ] ;
    prob.barc.subl=[prob.barc.subl J' ];
    prob.barc.val=[prob.barc.val V' ];
    clear I J V
    
% The barA matrices are formed from the rows of A
    for i=1:m
        barA= mat(A(i,Kindex(j):(Kindex(j+1)-1)),K.s(j));% j references the variable number (K.s(j)) and i references the ith row of A
        [I,J,V] = find(tril(barA+barA')/2);
        prob.bara.subi=[prob.bara.subi i*ones(1,length(I))];
        prob.bara.subj=[prob.bara.subj j*ones(1,length(I))];
        prob.bara.subk=[prob.bara.subk I' ];
        prob.bara.subl=[prob.bara.subl J' ];
        prob.bara.val=[prob.bara.val V' ];
        clear I J V
    end
end

prob.a=sparse(lpA);
prob.c=sparse(lpc);
prob.bardim=K.s;

prob.blc=b; %
prob.buc=b;
prob.blx=sparse([-inf*ones(nf,1);zeros(nl,1)]);
prob.bux=[];


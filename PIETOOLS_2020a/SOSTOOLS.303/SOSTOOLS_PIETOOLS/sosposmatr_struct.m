function [sos,P] = sosposmatr_struct(sos,n,sp_pat,wscoeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatrixTools v.03 - sosposmatr_struct(sos,n,wscoeff)
% This program declares a symbolic positive scalar semidefinite matrix P of size nxn 
%
% INPUTS: 
% sos - The SOSprogram to which to append the variable
% n - size of matrix variable
% wscoeff - if wscoeff='wscoeff', elements of the matrix will be declared
% as pvars in the workspace
%
% NOTES:
% Distributed with DelayTOOLS
% Compatable with MULTIPOLY and SOSTOOLS as of June 2013
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MatrixTools - sosposmatr
% Release Notes:
% v.03 - deleted the requirement to pass in a random symbolic
% variable
%
% Other changes with this release are minor.
%
% version .03   M. Peet
% MMP 6/2/2013 - restructured to reduce overhead. Eliminating the need for 
% Z by direct assignment of variables without call to sossosvar. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Copyright (C)2013  M. Peet
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


% So for now we are just going to take the easy way - declare a large
% positive semidefinite matrix and only return the requested variables.
% This assumes that sossolve doesn't really need to know which variables
% are 

%%
% First add a positive matrix of the appropriate dimension. This is all
% repeated from sosposmatr. The only difference is that we don't actually
% construct the matrix, which is time-consuming.

if ~all(size(sp_pat)==n)
    error('sparsity pattern must be nxn')
end
if ~all(all(triu(sp_pat)==tril(sp_pat)'))
    error('sparsity pattern must be symmetric')
end
    
    % Add new variable
    sos.var.num = sos.var.num+1;
    var = sos.var.num;
    sos.var.type{var} = 'sos';
    sos.var.Z{var} = [];%makesparse(Z);
    sos.var.ZZ{var} = [];
    sos.var.T{var} = [];
    sos.var.idx{var+1} = sos.var.idx{var}+n^2;


    for i = 1:sos.expr.num
        sos.expr.At{i} = [sos.expr.At{i}; ...
            sparse(n^2,size(sos.expr.At{i},2))];
    end;
    sos.objective = [sos.objective; sparse(sos.var.idx{var+1}-sos.var.idx{var},1)];        % 01/07/02

%%   
% Now that the positive matrix has been added, we would like to only
% declare the relevant decision variables. The problem is that the size of
% the decision variable table should probably match the length of At(,2).
% Now sossolve doesn't use the size of the decvartable. Neither does
% sosgetsol. 

[lind1] = find(triu(sp_pat));

% This section of code is for if you DO want to store the names of
% decision variables which do not appear in the program. If so, uncomment
% the following three lines and comment out the three lines after
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    sos.decvartable = [sos.decvartable; cellstr([repmat('coeff_',n^2,1), int2str(sos.var.idx{var}-sos.var.idx{1}+[1:n^2]')])];
%int2str(sos.var.idx{var}-sos.var.idx{1}+lind1)
%[repmat('coeff_',length(lind1),1), int2str(sos.var.idx{var}-sos.var.idx{1}+lind1)];

% This section of code is for if you do not want to store the names of
% decision variables which do not appear in the program. However, to be
% safe, we do not eliminate these variables, only zero out there name in
% the decvartable. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp=cellstr(char(ones(n^2,1)));
temp(lind1)=cellstr([repmat('coeff_',length(lind1),1), strjust(int2str(sos.var.idx{var}-sos.var.idx{1}+lind1),'left')]);
sos.decvartable = [sos.decvartable; temp];

        
    
% I should probably optimize this too. But really, who uses this option?
    if nargin > 3 & wscoeff == 'wscoeff'
        for i = lind1%sos.var.idx{var}:sos.var.idx{var+1}-1
            pvar([temp{i}]);
            assignin('base',[temp{i}],eval([temp{i}]));
        end;
	end;
 
nvars=length(lind1);

%begin_n=sos.var.idx{length(sos.var.idx)-1};
nTT=n;

% Lets try and create this directly
% Lets eliminate uneccesary variables by enforcing the symmetry constraint.
%
matdim=[nTT,nTT];
degmat=speye(nvars);

varname = temp(lind1);
[Ilist, Jlist] = ind2sub([n n], lind1);
lind2=sub2ind([n n], Jlist, Ilist);
nterms=length(lind1);

coeff = spones(sparse(repmat([1:nterms]',2,1),[lind1;lind2],1,nterms,(n)^2));

P=polynomial(coeff,degmat,varname,matdim);


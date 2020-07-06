function [sos,P] = sosposmatr(sos,n,wscoeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MatrixTools v.03 - sosposmatr(sos,n,wscoeff)
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

% ZZ and Z are not actually used, so this is probably wasted time. Will
% delete later after careful check...

% pvar(sos.vartable{1});
%     ZSym=monomials(eval(sos.vartable{1}),[1:n]);
% 
% 
%     if isnumeric(ZSym) & ZSym == 1
% 	pvar ZSym;
%         ZSym = ZSym*0+1;
%     end;
%     
%     [dum,idx1,idx2] = intersect(ZSym.varname,sos.vartable);
%     Z = sparse(n,length(sos.vartable));
%     Z(:,idx2) = sparse(ZSym.degmat(:,idx1));
    lenZ = n;  

    % Error handling needs to be added here, e.g. if Z is not a
    % valid monomial vector.

    
    
    
    % Add new variable
    sos.var.num = sos.var.num+1;
    var = sos.var.num;
    sos.var.type{var} = 'sos';
    sos.var.Z{var} = [];%makesparse(Z);
%    [T,ZZ] = getconstraint(Z);    
    sos.var.ZZ{var} = [];
    sos.var.T{var} = [];
    sos.var.idx{var+1} = sos.var.idx{var}+n^2;


    for i = 1:sos.expr.num
        sos.expr.At{i} = [sos.expr.At{i}; ...
            sparse(n^2,size(sos.expr.At{i},2))];
    end;


% pvar(sos.vartable{1});
%     ZSym=monomials(eval(sos.vartable{1}),[1:n]);
% 
% 
%     if isnumeric(ZSym) & ZSym == 1
% 	pvar ZSym;
%         ZSym = ZSym*0+1;
%     end;
%     
%     [dum,idx1,idx2] = intersect(ZSym.varname,sos.vartable);
%     Z = sparse(n,length(sos.vartable));
%     Z(:,idx2) = sparse(ZSym.degmat(:,idx1));
%     lenZ = n;  
% 
%     % Error handling needs to be added here, e.g. if Z is not a
%     % valid monomial vector.
% 
%     
%     
%     
%     % Add new variable
%     sos.var.num = sos.var.num+1;
%     var = sos.var.num;
%     sos.var.type{var} = 'sos';
%     sos.var.Z{var} = makesparse(Z);
%     [T,ZZ] = getconstraint(Z);    
%     sos.var.ZZ{var} = ZZ;
%     sos.var.T{var} = T';
%     sos.var.idx{var+1} = sos.var.idx{var}+size(Z,1)^2;
% 
    %T
%lenZ
%size(sos.var.T{var},1)
    % Modify existing equations
%     % Modify existing equations
%     for i = 1:sos.expr.num
%         sos.expr.At{i} = [sos.expr.At{i}; ...
%             sparse(size(sos.var.T{var},1),size(sos.expr.At{i},2))];
%     end;
    % Modify existing objective
    sos.objective = [sos.objective; sparse(sos.var.idx{var+1}-sos.var.idx{var},1)];        % 01/07/02

    % Modify decision variable table
%    oldlen = length(sos.decvartable);

[lind1] = find(triu(ones(lenZ)));

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
temp=cellstr(char(ones(lenZ^2,1)));
temp(lind1)=cellstr([repmat('coeff_',length(lind1),1), int2str(sos.var.idx{var}-sos.var.idx{1}+lind1)]);
sos.decvartable = [sos.decvartable; temp];

%     
%     sos.decvartable = [sos.decvartable; cellstr(char(ones(lenZ^2,2)))];
%     [lind1] = find(triu(ones(lenZ)));
% 
%     for i = lind1'
%         sos.decvartable(oldlen+i) = {['coeff_',int2str(sos.var.idx{var}-sos.var.idx{1}+i)]};
%     end;    
        
    
% I should probably optimize this two. But really, who uses this option?
    if nargin > 2 & wscoeff == 'wscoeff'
        var = sos.var.num;
        for i = sos.var.idx{var}:sos.var.idx{var+1}-1
            pvar(['coeff_',int2str(i-sos.var.idx{1}+1)]);
            assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],eval(['coeff_',int2str(i-sos.var.idx{1}+1)]));
        end;
	end;

    
%     for i = 1:lenZ^2
%         sos.decvartable(oldlen+i) = {['coeff_',int2str(sos.var.idx{var}-sos.var.idx{1}+i)]};
%     end;    
%  
%         
%     
%     if nargin > 2 & wscoeff == 'wscoeff'
%         var = sos.var.num;
%         for i = sos.var.idx{var}:sos.var.idx{var+1}-1
%             pvar(['coeff_',int2str(i-sos.var.idx{1}+1)]);
%             assignin('base',['coeff_',int2str(i-sos.var.idx{1}+1)],eval(['coeff_',int2str(i-sos.var.idx{1}+1)]));
%         end;
% 	end;

   
    

begin_n=sos.var.idx{length(sos.var.idx)-1};
nTT=n;

% Lets try and create this directly
% Lets eliminate uneccesary variables by enforcing the symmetry constraint.
%
matdim=[nTT,nTT];

degmat=speye((nTT^2/2+nTT/2));
%[lind1] = find(triu(ones(n)));
varname = sos.decvartable(begin_n-1+lind1);
[Ilist, Jlist] = ind2sub([n n], lind1);
lind2=sub2ind([n n], Jlist, Ilist);
nterms=length(lind1);

coeff = spones(sparse(repmat([1:nterms]',2,1),[lind1;lind2],1,nterms,(n)^2));

P=polynomial(coeff,degmat,varname,matdim);

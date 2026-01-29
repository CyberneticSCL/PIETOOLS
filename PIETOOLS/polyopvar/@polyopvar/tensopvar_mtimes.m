function AB = tensopvar_mtimes(A,B)
% AB = tensopvar_mtimes(A,B) returns the 'tensopvar' object C representing
% coefficients of the product A*B. Assumes A,B are represented with a common basis.
%
% INPUTS
% - A:     polyopvar object; 
% - B:     polyopvar object;
%
% OUTPUS
% - AB:    tensopvar object representing coefficients of A*B.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - tensopvar_mtimes
%
% Copyright (C) 2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% CR, 01/28/2026: Initial coding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for empty cells in A.C.ops and B.C.ops, then replace by zero nopvars.
% Assumes A,B have a common basis; i.e., nz=numel(A.C.ops)=numel(B.C.ops)
% and degmat=A.degmat=B.degmat.
nz=numel(A.C.ops); % 1xnz
degmat=A.degmat;

for i=1:nz
    Ai=A.C.ops{i}; % mi x di 
    Bi=B.C.ops{i}; % mi x di 
    [mi_a,~]=size(Ai); % [0,0] if Ai is empty.
    [mi_b,~]=size(Bi); % [0,0] if Bi is empty.
    if mi_a==0
        di=sum(A.degmat(i,:)); % tells us how many nopvar coeffs to add.
        C=zero_nopvar(A);
        A.C.ops{i}=repmat({C},1,di);
    end
    if mi_b==0
        di=sum(B.degmat(i,:)); % tells us how many nopvar coeffs to add.
        C=zero_nopvar(B);
        B.C.ops{i}=repmat({C},1,di);
    end
end

% Iterate through C{1:nz^2}{1:mi*qj,1:di+dj} and assign corresponding operators 
% from A and B in the required order.
C=cell(1,nz^2);
for i=1:nz
    Ai=A.C.ops{i}; % mi x di 
    [mi,di]=size(Ai);
    for j=1:nz
        Bj=B.C.ops{j}; % qj x di
        [qj,dj]=size(Bj);
        for k=1:mi
            for l=1:qj
                alpha=(i-1)*nz+j;
                beta=(k-1)*qj+l;
                C{alpha}=cell(mi*qj,di+dj); % what we expect it to equal.
                p=1; % starting column index of degmat.
                count_a=1; count_b=1; count_c=1; % track operator positions.
                while p<size(degmat,2)
                    if degmat(i,p)==0 && degmat(j,p)==0 % move to next variable.
                        p=p+1; 
                    elseif degmat(i,p)==0 && degmat(j,p)>0 % update C with B, then move to next variable.
                        deg=degmat(j,p);
                        [C{alpha}{beta,count_c:count_c+deg-1}]=Bj{l,count_b:count_b+deg-1};
                        count_c=count_c+deg;
                        count_b=count_b+deg;
                        p=p+1;
                    elseif degmat(i,p)>0 && degmat(j,p)==0 % update C with A, then move to next variable.
                        deg=degmat(i,p);
                        [C{alpha}{beta,count_c:count_c+deg-1}]=Ai{k,count_a:count_a+deg-1};
                        count_c=count_c+deg;
                        count_a=count_a+deg;
                        p=p+1;
                    else % update C with A and B, then move to next variable.
                        deg=degmat(i,p);
                        [C{alpha}{beta,count_c:count_c+deg-1}]=Ai{k,count_a:count_a+deg-1};
                        count_c=count_c+deg;
                        count_a=count_a+deg;
                        deg=degmat(j,p);
                        [C{alpha}{beta,count_c:count_c+deg-1}]=Bj{l,count_b:count_b+deg-1};
                        count_c=count_c+deg;
                        count_b=count_b+deg;
                        p=p+1;
                    end
                end
            end
        end
    end
end

AB = tensopvar();
AB.ops = C;

end
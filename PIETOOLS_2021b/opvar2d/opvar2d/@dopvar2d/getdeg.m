function [T] = getdeg(P,fprint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [T] = getdeg(P,fprint) provides the min and max degree of the monomials present in
% the components of a PI operator. The result returned is a 2D array with
% entries described as below.
% For a PI operator with components P,Q1(s),Q2(s),R0(s),R1(s,t),R2(s,t), 
% the function returns
% the following table.
% T = [  min_P_s  minP_t  max_P_s  max_P_t;
%        min_Q1_s minQ1_t max_Q1_s max_Q1_t;
%        min_Q2_s minQ2_t max_Q2_s max_Q2_t;
%        min_R0_s minR0_t max_R0_s max_R0_t;
%        min_R1_s minR1_t max_R1_s max_R1_t;
%        min_R2_s minR2_t max_R2_s max_R2_t]
% where the variables naming follows boundtype_componentname_variable.
% 
% INPUT
% P : opvar class object
% fprint: flag if set to 1 the degree bounds will be printed
% 
% OUTPUT
% T : Array containing min and max degree of monomials in the components of P
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - getdeg
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 7_26_2019
% Changed datatype of var1, var2 to include pvars and chars - SS 9/6/2019
% Updated for dopvar implementation - DJ 12/27/2021

if nargin==1
    fprint=0;
end

if ~isa(P,'dopvar')
    error('Input must be a dopvar variable.');
end
fset = {'P','Q1','Q2'};
k=1;
T= zeros(6,4);

if isa(P.var1,'polynomial')
    var1 = P.var1.varname;
elseif isa(P.var1,'char')
    var1 = P.var1;
end
if isa(P.var2,'polynomial')
    var2 = P.var2.varname;
elseif isa(P.var2,'char')
    var2 = P.var2;
end


for i=fset
if isa(P.(i{:}),'polynomial') || isa(P.(i{:}),'dpvar')
    idx1= find(strcmp(P.(i{:}).varname,var1));
    idx2= find(strcmp(P.(i{:}).varname,var2));
    deg = P.(i{:}).degmat(:,[idx1,idx2]);
    if ~isempty(idx1)
        T(k,[1,3]) = [min(deg(:,1)) max(deg(:,1))];
        if ~isempty(idx2)
            T(k,[2,4]) = [min(deg(:,2)) max(deg(:,2))];
        end
    else
        if ~isempty(idx2)
            T(k,[2,4]) = [min(deg(:,1)) max(deg(:,1))];
        end
    end
end
k=k+1;
end

fset = {'R0','R1','R2'};
for i=fset
if isa(P.R.(i{:}),'polynomial') || isa(P.R.(i{:}),'dpvar')
    idx1= find(strcmp(P.R.(i{:}).varname,var1));
    idx2= find(strcmp(P.R.(i{:}).varname,var2));
    deg = P.R.(i{:}).degmat(:,[idx1,idx2]);
    if ~isempty(idx1)
        T(k,[1,3]) = [min(deg(:,1)) max(deg(:,1))];
        if ~isempty(idx2)
            T(k,[2,4]) = [min(deg(:,2)) max(deg(:,2))];
        end
    else
        if ~isempty(idx2)
            T(k,[2,4]) = [min(deg(:,1)) max(deg(:,1))];
        end
    end
end
k=k+1;
end
T = full(T);
if fprint==1
min_var1 = T(:,1);
min_var2 = T(:,2);
max_var1 = T(:,3);
max_var2 = T(:,4);
Components = {'P','Q1','Q2','R.R0','R.R1','R.R2'}.';
table(Components,min_var1,min_var2,max_var1,max_var2)
end
end

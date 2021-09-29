function [T] = getdeg(P,fprint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [T] = getdeg(P,fprint) provides the min and max degree of the monomials present in
% the components of a decision PI operator. The result returned is a 2D array with
% entries described as below.
% For a PI operator with components 
%   R00 , R0x , R0y , R02;
%   Rx0 , Rxx , Rxy , Rx2;
%   Ry0 , Ryx , Ryy , Ry2;
%   R20 , R2x , R2y , R22;
% the function returns
% the following table.
% T = [ min_R00_ss1         minR00_tt1      min_R00_ss2         min_R00_tt2;
%       min_R0x_ss1         minR0x_tt1      min_R0x_ss2         min_R0x_tt2;
%       min_R0y_ss1         minR0y_tt1      min_R0y_ss2         min_R0y_tt2;
%       min_R02_ss1         minR02_tt1      min_R02_ss2         min_R02_tt2;
%       min_Rx0_ss1         minRx0_tt1      min_Rx0_ss2         min_Rx0_tt2;
%       min_Rxx{1}_ss1      minRxx{1}_tt1   min_Rxx{1}_ss2      min_Rxx{1}_tt2;
%       min_Rxx{2}_ss1      minRxx{2}_tt1   min_Rxx{2}_ss2      min_Rxx{2}_tt2;
%       min_Rxx{3}_ss1      minRxx{3}_tt1   min_Rxx{3}_ss2      min_Rxx{3}_tt2;
%       min_Rxy_ss1         minRxy_tt1      min_Rxy_ss2         min_Rxy_tt2;
%       min_Rx2{1}_ss1      minRx2{1}_tt1   min_Rx2{1}_ss2      min_Rx2{1}_tt2;
%       min_Rx2{2}_ss1      minRx2{2}_tt1   min_Rx2{2}_ss2      min_Rx2{2}_tt2;
%       min_Rx2{3}_ss1      minRx2{3}_tt1   min_Rx2{3}_ss2      min_Rx2{3}_tt2;
%       min_Ry0_ss1         minRy0_tt1      min_Ry0_ss2         min_Ry0_tt2;
%       min_Ryx_ss1         minRyx_tt1      min_Ryx_ss2         min_Ryx_tt2;
%       min_Ryy{1}_ss1      minRyy{1}_tt1   min_Ryy{1}_ss2      min_Ryy{1}_tt2;
%       min_Ryy{2}_ss1      minRyy{2}_tt1   min_Ryy{2}_ss2      min_Ryy{2}_tt2;
%       min_Ryy{3}_ss1      minRyy{3}_tt1   min_Ryy{3}_ss2      min_Ryy{3}_tt2;
%       min_Ry2{1}_ss1      minRy2{1}_tt1   min_Ry2{1}_ss2      min_Ry2{1}_tt2;
%       min_Ry2{2}_ss1      minRy2{2}_tt1   min_Ry2{2}_ss2      min_Ry2{2}_tt2;
%       min_Ry2{3}_ss1      minRy2{3}_tt1   min_Ry2{3}_ss2      min_Ry2{3}_tt2;
%       min_R20_ss1         minR20_tt1      min_R20_ss2         min_R20_tt2;
%       min_R2x{1}_ss1      minR2x{1}_tt1   min_R2x{1}_ss2      min_R2x{1}_tt2;
%       min_R2x{2}_ss1      minR2x{2}_tt1   min_R2x{2}_ss2      min_R2x{2}_tt2;
%       min_R2x{3}_ss1      minR2x{3}_tt1   min_R2x{3}_ss2      min_R2x{3}_tt2;
%       min_R2y{1}_ss1      minR2y{1}_tt1   min_R2y{1}_ss2      min_R2y{1}_tt2;
%       min_R2y{2}_ss1      minR2y{2}_tt1   min_R2y{2}_ss2      min_R2y{2}_tt2;
%       min_R2y{3}_ss1      minR2y{3}_tt1   min_R2y{3}_ss2      min_R2y{3}_tt2;
%       min_R22{1,1}_ss1    minR22{1,1}_tt1 min_R22{1,1}_ss2    min_R22{1,1}_tt2;
%       min_R22{2,1}_ss1    minR22{2,1}_tt1 min_R22{2,1}_ss2    min_R22{2,1}_tt2;
%       min_R22{3,1}_ss1    minR22{3,1}_tt1 min_R22{3,1}_ss2    min_R22{3,1}_tt2;
%       min_R22{1,2}_ss1    minR22{1,2}_tt1 min_R22{1,2}_ss2    min_R22{1,2}_tt2;
%       min_R22{1,3}_ss1    minR22{1,3}_tt1 min_R22{1,3}_ss2    min_R22{1,3}_tt2;
%       min_R22{2,2}_ss1    minR22{2,2}_tt1 min_R22{2,2}_ss2    min_R22{2,2}_tt2;
%       min_R22{3,2}_ss1    minR22{3,2}_tt1 min_R22{3,2}_ss2    min_R22{3,2}_tt2;
%       min_R22{2,3}_ss1    minR22{2,3}_tt1 min_R22{2,3}_ss2    min_R22{2,3}_tt2;
%       min_R22{3,3}_ss1    minR22{3,3}_tt1 min_R22{3,3}_ss2    min_R22{3,3}_tt2]
%
% 
% INPUT
% P : dopvar2d class object
% fprint: flag if set to 1 the degree bounds will be printed
% 
% OUTPUT
% T : Array containing min and max degree of monomials in the components of P
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - getdeg
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 12_07_2021

if nargin==1
    fprint=0;
end

if ~isa(P,'dopvar2d')
    error('Input must be a dopvar2d variable.');
end
fset = {'R00','R0x','R0y','R02','Rx0','Rxy','Ry0','Ryx','R20'};
kset = {1,2,3,4,5,9,13,14,21};
fsetx = {'Rxx','Rx2','R2x'};
ksetx = {(6:8),(10:12),(22:24)};
fsety = {'Ryy','Ry2','R2y'};
ksety = {(15:17),(18:20),(25:27)};
fset2 = {'R22'};
kset2 = {[28,31,32;29,33,35;30,34,36]};
k=1;
ftot = length(fset)+3*length(fsetx)+3*length(fsety)+9*length(fset2);
T = zeros(ftot,8);

% ss1
if isa(P.var1(1),'polynomial')
    var11 = P.var1(1).varname;
elseif isa(P.var1,'char')
    var11 = P.var1(1);
end
% tt1
if isa(P.var2(1),'polynomial')
    var12 = P.var2(1).varname;
elseif isa(P.var2,'char')
    var12 = P.var2(1);
end
% ss2
if isa(P.var1(2),'polynomial')
    var21 = P.var1(2).varname;
elseif isa(P.var1,'char')
    var21 = P.var1(2);
end
% tt2
if isa(P.var2(2),'polynomial')
    var22 = P.var2(2).varname;
elseif isa(P.var2,'char')
    var22 = P.var2(2);
end


for f=1:length(fset)
k = kset{f};
    if isa(P.(fset{f}),'polynomial') || isa(P.(fset{f}),'dpvar')
        idx1 = find(strcmp(P.(fset{f}).varname,var11));
        idx2 = find(strcmp(P.(fset{f}).varname,var12));
        idx3 = find(strcmp(P.(fset{f}).varname,var21));
        idx4 = find(strcmp(P.(fset{f}).varname,var22));
        idx = [idx1,idx2,idx3,idx4];
        deg = P.(fset{f}).degmat(:,[idx1,idx2,idx3,idx4]);

        if ~isempty(idx1)
            i1 = find(idx==idx1);
            T(k,[1,5]) = [min(deg(:,i1)) max(deg(:,i1))];
        end
        if ~isempty(idx2)
            i1 = find(idx==idx2);
            T(k,[2,6]) = [min(deg(:,i1)) max(deg(:,i1))];
        end
        if ~isempty(idx3)
            i1 = find(idx==idx3);
            T(k,[3,7]) = [min(deg(:,i1)) max(deg(:,i1))];
        end
        if ~isempty(idx4)
            i1 = find(idx==idx4);
            T(k,[4,8]) = [min(deg(:,i1)) max(deg(:,i1))];
        end  
    end
end

fsetxy = [fsetx, fsety];
ksetxy = [ksetx, ksety];
for f=1:length(fsetxy)
    for i=1:3
        k = ksetxy{f}(i);
        if isa(P.(fsetxy{f}){i},'polynomial') || isa(P.(fsetxy{f}){i},'dpvar')
            idx1 = find(strcmp(P.(fsetxy{f}){i}.varname,var11));
            idx2 = find(strcmp(P.(fsetxy{f}){i}.varname,var12));
            idx3 = find(strcmp(P.(fsetxy{f}){i}.varname,var21));
            idx4 = find(strcmp(P.(fsetxy{f}){i}.varname,var22));
            idx = [idx1,idx2,idx3,idx4];
            deg = P.(fsetxy{f}){i}.degmat(:,[idx1,idx2,idx3,idx4]);

            if ~isempty(idx1)
                i1 = find(idx==idx1);
                T(k,[1,5]) = [min(deg(:,i1)) max(deg(:,i1))];
            end
            if ~isempty(idx2)
                i1 = find(idx==idx2);
                T(k,[2,6]) = [min(deg(:,i1)) max(deg(:,i1))];
            end
            if ~isempty(idx3)
                i1 = find(idx==idx3);
                T(k,[3,7]) = [min(deg(:,i1)) max(deg(:,i1))];
            end
            if ~isempty(idx4)
                i1 = find(idx==idx4);
                T(k,[4,8]) = [min(deg(:,i1)) max(deg(:,i1))];
            end  
        end
    end
end

for f=1:length(fset2)
    for i=1:3
    for j=1:3
        k = kset2{f}(i,j);
        if isa(P.(fset2{f}){i,j},'polynomial') || isa(P.(fset2{f}){i,j},'dpvar')
            idx1 = find(strcmp(P.(fset2{f}){i,j}.varname,var11));
            idx2 = find(strcmp(P.(fset2{f}){i,j}.varname,var12));
            idx3 = find(strcmp(P.(fset2{f}){i,j}.varname,var21));
            idx4 = find(strcmp(P.(fset2{f}){i,j}.varname,var22));
            idx = [idx1,idx2,idx3,idx4];
            deg = P.(fset2{f}){i,j}.degmat(:,[idx1,idx2,idx3,idx4]);

            if ~isempty(idx1)
                i1 = find(idx==idx1);
                T(k,[1,5]) = [min(deg(:,i1)) max(deg(:,i1))];
            end
            if ~isempty(idx2)
                i1 = find(idx==idx2);
                T(k,[2,6]) = [min(deg(:,i1)) max(deg(:,i1))];
            end
            if ~isempty(idx3)
                i1 = find(idx==idx3);
                T(k,[3,7]) = [min(deg(:,i1)) max(deg(:,i1))];
            end
            if ~isempty(idx4)
                i1 = find(idx==idx4);
                T(k,[4,8]) = [min(deg(:,i1)) max(deg(:,i1))];
            end  
        end
    end
    end
end



T = full(T);
if fprint==1
    fsetfull = {'R00','R0x','R0y','R02',...
                'Rx0','Rxx{1}','Rxx{2}','Rxx{3}','Rxy','Rx2{1}','Rx2{2}','Rx2{3}',...
                'Ry0','Ryx','Ryy{1}','Ryy{2}','Ryy{3}','Ry2{1}','Ry2{2}','Ry2{3}',...
                'R20','R2x{1}','R2x{2}','R2x{3}','R2y{1}','R2y{2}','R2y{3}',...
                'R22{1,1}','R22{2,1}','R22{3,1}','R22{1,2}','R22{1,3}',...
                'R22{2,2}','R22{3,2}','R22{2,3}','R22{3,3}'}';
    min_var11 = T(:,1);
    min_var12 = T(:,2);
    min_var21 = T(:,3);
    min_var22 = T(:,4);
    max_var11 = T(:,5);
    max_var12 = T(:,6);
    max_var21 = T(:,7);
    max_var22 = T(:,8);
    table(fsetfull,min_var11,min_var12,min_var21,min_var22,max_var11,max_var12,max_var21,max_var22)
end

end
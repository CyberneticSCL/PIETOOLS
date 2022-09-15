function display(T,name)
% display(T,name) prints the components of the opvar object P
% 
% INPUTS:
% P : An opvar variable
% name: name of variable used 
% 
% OUTPUTS:
% command line display of elements of opvar
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS  - 6_20_2020


slimit = [8,8];

Pbig = 0; Q1big = 0; Q2big = 0; R0big = 0; R1big = 0; R2big = 0;

opts = {'P','Q2','R0','R1','R2'};


if ~isa(T,'dopvar')
    error('Input must be an opvar variable.');
end


P = T.P; Q1 = T.Q1; Q2 = T.Q2; R0 = T.R.R0; R1 = T.R.R1; R2 = T.R.R2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the matrix dimensions of the components are too large to display
if any(size(P)>slimit)
    Pbig = 1;
end
if any(size(Q1)>slimit)
    Q1big = 1;
end
if any(size(Q2)>slimit)
    Q2big = 1;
end
if any(size(R0)>slimit)
    R0big = 1;
end
if any(size(R1)>slimit)
    R1big = 1;
end
if any(size(R2)>slimit)
    R2big = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lets get output string together
[Pstr, Pbig] = getstringout(P,Pbig);
[Q1str, Q1big] = getstringout(Q1,Q1big);
[Q2str, Q2big] = getstringout(Q2,Q2big);
[R0str, R0big] = getstringout(R0,R0big);
[R1str, R1big] = getstringout(R1,R1big);
[R2str, R2big] = getstringout(R2,R2big);




if Pbig
    Pstr = "Too big to display. Use ans.P";
end
if Q1big
    Q1str = "Too big to display. Use ans.Q1";
end
if Q2big
    Q2str = "Too big to display. Use ans.Q2";
end
if R0big
    R0str = "Too big to display. Use ans.R.R0";
end
if R1big
    R1str = "Too big to display. Use ans.R.R1";
end
if R2big
    R2str = "Too big to display. Use ans.R.R2";
end

% display RL2 to RL2 part
disp([name,'=']);
fprintf("\n");
Rmx = max([size(Pstr,2),size(Q1str,2)]);
if size(Pstr,2)<Rmx
    Pstr(end+1:Rmx) = " ";
end
if size(Q1str,2)<Rmx
    Q1str(end+1:Rmx) = " ";
end

Plen = max(arrayfun(@(x) strlength(x),Pstr));
Qlen = max(arrayfun(@(x) strlength(x),Q2str));
for i = 1:Rmx
    fprintf("%5s %*s | %s \n"," ", max(Plen,Qlen),Pstr(i),Q1str(i));
end
PQ = Pstr+"|"+Q1str; 
TRstr = name+".R"; TRstr(end+1:size(Q2str,2)) = " ";
QR = Q2str+"|"+TRstr;
divlen = max(max(arrayfun(@(x) strlength(x),PQ)),max(arrayfun(@(x) strlength(x),QR))); 
div = "";
for i=1:divlen+3
div = div+"-";
end
fprintf("%5s %s\n%"," ",div);
for i = 1:size(Q2str,2)
    fprintf("%5s %*s | %s \n"," ",max(Plen,Qlen),Q2str(i),TRstr(i));
end


fprintf("\n");
% disp L2 to L2 part
disp([name,'.R=']);
fprintf("\n");
Rmx = max([size(R0str,2),size(R1str,2),size(R2str,2)]);
if size(R0str,2)<Rmx
    R0str(end+1:Rmx) = " ";
end
if size(R1str,2)<Rmx
    R1str(end+1:Rmx) = " ";
end
if size(R2str,2)<Rmx
    R2str(end+1:Rmx) = " ";
end
R0rowmax = max(arrayfun(@(x) strlength(x),R0str));
R1rowmax = max(arrayfun(@(x) strlength(x),R1str));
R2rowmax = max(arrayfun(@(x) strlength(x),R2str));
for i = 1:Rmx
    fprintf("    %*s | %*s | %*s \n", R0rowmax, R0str(i),R1rowmax, R1str(i),R2rowmax,R2str(i));
end

% disp([name,'.dim=']);
% T.dim
% disp([name,'.I=']);
% T.I

end
function [str,big] = getstringout(prop,big)
llimit = 80; 
if isempty(prop)
str = "[]"; 
else
if ~big
    for i=1:size(prop,1)
        row ="";
        for j=1:size(prop,2)
            try 
                temp = evalc('prop(i,j)');
            catch
                str="";
                big = 1;
                return
            end
            temp = convertCharsToStrings(temp);
            temp = regexprep(temp,'[\n\r]','');
            temp = regexprep(temp,' ','');
            temp = regexprep(temp,'ans=','');
            row = row+temp+",";
        end
        row = convertStringsToChars(row);
        row = row(1:end-1);
        row = convertCharsToStrings(row);
        if strlength(row)<=llimit
            str(i) = "["+row+"]";
        else
            str="";
            big = 1;
        end
    end
else
    str="";
end
end
end
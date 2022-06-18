function display(T,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display(T,name) prints the components of the dopvar2d object P
% 
% INPUTS:
% P : A dopvar2d variable
% name: name of variable used 
% 
% OUTPUTS:
% command line display of elements of opvar
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07_12_2021  


slimit = [8,8];

if ~isa(T,'dopvar2d')
    error('Input must be a dopvar2d variable.');
end

% Collect the functions
R00 = T.R00;    R0x = T.R0x;    R0y = T.R0y;    R02 = T.R02;
Rx0 = T.Rx0;    Rxx = T.Rxx;    Rxy = T.Rxy;    Rx2 = T.Rx2;
Ry0 = T.Ry0;    Ryx = T.Ryx;    Ryy = T.Ryy;    Ry2 = T.Ry2;
R20 = T.R20;    R2x = T.R2x;    R2y = T.R2y;    R22 = T.R22;


R00big = 0; R0xbig = 0;             R0ybig = 0;             R02big = 0; 
Rx0big = 0; Rxxbig = zeros(3,1);    Rxybig = 0;             Rx2big = zeros(3,1);
Ry0big = 0; Ryxbig = 0;             Ryybig = zeros(1,3);    Ry2big = zeros(1,3);
R20big = 0; R2xbig = zeros(3,1);    R2ybig = zeros(1,3);    R22big = zeros(3,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the matrix dimensions of the components are too large to display
if any(size(R00)>slimit)
    R00big = 1;
end
if any(size(R0x)>slimit)
    R0xbig = 1;
end
if any(size(R0y)>slimit)
    R0ybig = 1;
end
if any(size(R02)>slimit)
    R02big = 1;
end
%
if any(size(Rx0)>slimit)
    Rx0big = 1;
end
if any(size(Rxy)>slimit)
    Rxybig = 1;
end
%
if any(size(Ry0)>slimit)
    Ry0big = 1;
end
if any(size(Ryx)>slimit)
    Ryxbig = 1;
end
%
if any(size(R20)>slimit)
    R20big = 1;
end

for i=1:3
    if any(size(Rxx{i,1})>slimit)
        Rxxbig(i,1) = 1;
    end
    if any(size(Rx2{i,1})>slimit)
        Rx2big(i,1) = 1;
    end
    if any(size(Ryy{1,i})>slimit)
        Ryybig(1,i) = 1;
    end
    if any(size(Ry2{1,i})>slimit)
        Ry2big(1,i) = 1;
    end

    if any(size(R2x{i,1})>slimit)
        R2xbig(i,1) = 1;
    end
    if any(size(R2y{1,i})>slimit)
        R2ybig(1,i) = 1;
    end
    
    for j=1:3
        if any(size(R22{i,j})>slimit)
            R22big(i,j) = 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lets get output string together
[R00str, R00big] = getstringout(R00,R00big);
[R0xstr, R0xbig] = getstringout(R0x,R0xbig);
[R0ystr, R0ybig] = getstringout(R0y,R0ybig);
[R02str, R02big] = getstringout(R02,R02big);
%
[Rx0str, Rx0big] = getstringout(Rx0,Rx0big);
[Rxxstr_o, Rxxbigo] = getstringout(Rxx{1},Rxxbig(1));
[Rxxstr_a, Rxxbiga] = getstringout(Rxx{2},Rxxbig(2));
[Rxxstr_b, Rxxbigb] = getstringout(Rxx{3},Rxxbig(3));
Rxxbig = [Rxxbigo; Rxxbiga; Rxxbigb];
[Rxystr, Rxybig] = getstringout(Rxy,Rxybig);
[Rx2str_o, Rx2bigo] = getstringout(Rx2{1},Rx2big(1));
[Rx2str_a, Rx2biga] = getstringout(Rx2{2},Rx2big(2));
[Rx2str_b, Rx2bigb] = getstringout(Rx2{3},Rx2big(3));
Rx2big = [Rx2bigo; Rx2biga; Rx2bigb];
%
[Ry0str, Ry0big] = getstringout(Ry0,Ry0big);
[Ryxstr, Ryxbig] = getstringout(Ryx,Ryxbig);
[Ryystr_o, Ryybigo] = getstringout(Ryy{1},Ryybig(1));
[Ryystr_a, Ryybiga] = getstringout(Ryy{2},Ryybig(2));
[Ryystr_b, Ryybigb] = getstringout(Ryy{3},Ryybig(3));
Ryybig = [Ryybigo, Ryybiga, Ryybigb];
[Ry2str_o, Ry2bigo] = getstringout(Ry2{1},Ry2big(1));
[Ry2str_a, Ry2biga] = getstringout(Ry2{2},Ry2big(2));
[Ry2str_b, Ry2bigb] = getstringout(Ry2{3},Ry2big(3));
Ry2big = [Ry2bigo, Ry2biga, Ry2bigb];
%
[R20str, R20big] = getstringout(R20,R20big);
[R2xstr_o, R2xbigo] = getstringout(R2x{1},R2xbig(1));
[R2xstr_a, R2xbiga] = getstringout(R2x{2},R2xbig(2));
[R2xstr_b, R2xbigb] = getstringout(R2x{3},R2xbig(3));
R2xbig = [R2xbigo; R2xbiga; R2xbigb];
[R2ystr_o, R2ybigo] = getstringout(R2y{1},R2ybig(1));
[R2ystr_a, R2ybiga] = getstringout(R2y{2},R2ybig(2));
[R2ystr_b, R2ybigb] = getstringout(R2y{3},R2ybig(3));
R2ybig = [R2ybigo, R2ybiga, R2ybigb];
R22str = cell(3,3);
for i=1:3
    for j=1:3
        [Rstr_ij, Rbig_ij] = getstringout(R22{i,j},R22big(i,j));
        R22str{i,j} = Rstr_ij;
        R22big(i,j) = Rbig_ij;
        if Rbig_ij
            R22str{i,j} = "Too big to display. Use T.R22{"+num2str(i)+","+num2str(j)+"}";
        end
    end
end


if R00big
    R00str = "Too big to display. Use T.R00";
end
if R0xbig
    R0xstr = "Too big to display. Use T.R0x";
end
if R0ybig
    R0ystr = "Too big to display. Use T.R0y";
end
if R02big
    R02str = "Too big to display. Use T.R02";
end
%
if Rx0big
    Rx0str = "Too big to display. Use T.Rx0";
end
if Rxxbig(1)
    Rxxstr_o = "Too big to display. Use T.Rxx{1,1}";
end
if Rxxbig(2)
    Rxxstr_a = "Too big to display. Use T.Rxx{2,1}";
end
if Rxxbig(3)
    Rxxstr_b = "Too big to display. Use T.Rxx{3,1}";
end
if Rxybig
    Rxystr = "Too big to display. Use T.Rxy";
end
if Rx2big(1)
    Rx2str_o = "Too big to display. Use T.Rx2{1,1}";
end
if Rx2big(2)
    Rx2str_a = "Too big to display. Use T.Rx2{2,1}";
end
if Rx2big(3)
    Rx2str_b = "Too big to display. Use T.Rx2{3,1}";
end
%
if Ry0big
    Ry0str = "Too big to display. Use T.Ry0";
end
if Ryxbig
    Ryxstr = "Too big to display. Use T.Ryx";
end
if Ryybig(1)
    Ryystr_o = "Too big to display. Use T.Ryy{1,1}";
end
if Ryybig(2)
    Ryystr_a = "Too big to display. Use T.Ryy{1,2}";
end
if Ryybig(3)
    Ryystr_b = "Too big to display. Use T.Ryy{1,3}";
end
if Ry2big(1)
    Ry2str_o = "Too big to display. Use T.Ry2{1,1}";
end
if Ry2big(2)
    Ry2str_a = "Too big to display. Use T.Ry2{1,2}";
end
if Ry2big(3)
    Ry2str_b = "Too big to display. Use T.Ry2{1,3}";
end
%
if R20big
    R20str = "Too big to display. Use T.R20";
end
if R2xbig(1)
    R2xstr_o = "Too big to display. Use T.R2x{1,1}";
end
if R2xbig(2)
    R2xstr_a = "Too big to display. Use T.R2x{2,1}";
end
if R2xbig(3)
    R2xstr_b = "Too big to display. Use T.R2x{3,1}";
end
if R2ybig(1)
    R2ystr_o = "Too big to display. Use T.R2y{1,1}";
end
if R2ybig(2)
    R2ystr_a = "Too big to display. Use T.R2y{1,2}";
end
if R2ybig(3)
    R2ystr_b = "Too big to display. Use T.R2y{1,3}";
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Display non-subdivided parts
% First row
disp([name,' =']);
fprintf("\n");
Rmx0 = max([size(R00str,2),size(R0xstr,2),size(R0ystr,2),size(R02str,2)]);
if size(R00str,2)<Rmx0
    R00str(end+1:Rmx0) = " ";
end
if size(R0xstr,2)<Rmx0
    R0xstr(end+1:Rmx0) = " ";
end
if size(R0ystr,2)<Rmx0
    R0ystr(end+1:Rmx0) = " ";
end
if size(R02str,2)<Rmx0
    R02str(end+1:Rmx0) = " ";
end

Rmxx = max([size(Rx0str,2),size(Rxystr,2)]);
if size(Rx0str,2)<Rmxx
    Rx0str(end+1:Rmxx) = " ";
end
if size(Rxystr,2)<Rmxx
    Rxystr(end+1:Rmxx) = " ";
end

Rmxy = max([size(Ry0str,2),size(Ryxstr,2)]);
if size(Ry0str,2)<Rmxy
    Ry0str(end+1:Rmxy) = " ";
end
if size(Ryxstr,2)<Rmxy
    Ryxstr(end+1:Rmxy) = " ";
end

%R0 = R00str+"|"+R0xstr+"|"+R0ystr+"|"+R02str; 

TRxxstr = name+".Rxx";  TRxxstr(end+1:Rmxx) = " ";
TRx2str = name+".Rx2";  TRx2str(end+1:Rmxx) = " ";
%Rx = Rx0str+"|"+TRxxstr+"|"+Rxystr+"|"+TRx2str;

TRyystr = name+".Ryy";  TRyystr(end+1:Rmxy) = " ";
TRy2str = name+".Ry2";  TRy2str(end+1:Rmxy) = " ";
%Ry = Ry0str+"|"+Ryxstr+"|"+TRyystr+"|"+TRy2str;

TR2xstr = name+".R2x";  TR2xstr(end+1:size(R20str,2)) = " ";
TR2ystr = name+".R2y";  TR2ystr(end+1:size(R20str,2)) = " ";
TR22str = name+".R22";  TR22str(end+1:size(R20str,2)) = " ";
%R2 = R20str+"|"+TR2xstr+"|"+TR2ystr+"|"+TR22str;

R00len = max(arrayfun(@(x) strlength(x),R00str));
Rx0len = max(arrayfun(@(x) strlength(x),Rx0str));
Ry0len = max(arrayfun(@(x) strlength(x),Ry0str));
R20len = max(arrayfun(@(x) strlength(x),R20str));
mlen0 = max([R00len,Rx0len,Ry0len,R20len]);
R0xlen = max(arrayfun(@(x) strlength(x),R0xstr));
Rxxlen = max(arrayfun(@(x) strlength(x),TRxxstr));
Ryxlen = max(arrayfun(@(x) strlength(x),Ryxstr));
R2xlen = max(arrayfun(@(x) strlength(x),TR2xstr));
mlenx = max([R0xlen,Rxxlen,Ryxlen,R2xlen]);
R0ylen = max(arrayfun(@(x) strlength(x),R0ystr));
Rxylen = max(arrayfun(@(x) strlength(x),Rxystr));
Ryylen = max(arrayfun(@(x) strlength(x),TRyystr));
R2ylen = max(arrayfun(@(x) strlength(x),TR2ystr));
mleny = max([R0ylen,Rxylen,Ryylen,R2ylen]);
R02len = max(arrayfun(@(x) strlength(x),R02str));
Rx2len = max(arrayfun(@(x) strlength(x),TRx2str));
Ry2len = max(arrayfun(@(x) strlength(x),TRy2str));
R22len = max(arrayfun(@(x) strlength(x),TR22str));
mlen2 = max([R02len,Rx2len,Ry2len,R22len]);

divlen = mlen0 + mlenx + mleny + mlen2;
div = "";
for i=1:divlen+3*3
div = div+"-";
end

for i = 1:Rmx0
    fprintf("%3s %*s | %*s | %*s | %*s \n"," ", mlen0,R00str(i),mlenx,R0xstr(i),mleny,R0ystr(i),mlen2,R02str(i));
end
fprintf("%3s %s\n%"," ",div);
for i = 1:Rmxx
    fprintf("%3s %*s | %*s | %*s | %*s \n"," ",mlen0,Rx0str(i),mlenx,TRxxstr(i),mleny,Rxystr(i),mlen2,TRx2str(i));
end
fprintf("%3s %s\n%"," ",div);
for i = 1:Rmxy
    fprintf("%3s %*s | %*s | %*s | %*s \n"," ",mlen0,Ry0str(i),mlenx,Ryxstr(i),mleny,TRyystr(i),mlen2,TRy2str(i));
end
fprintf("%3s %s\n%"," ",div);
for i = 1:size(R20str,2)
    fprintf("%3s %*s | %*s | %*s | %*s \n"," ",mlen0,R20str(i),mlenx,TR2xstr(i),mleny,TR2ystr(i),mlen2,TR22str(i));
end


%%%%% Display subcomponents
%%% Rxx
fprintf("\n");
disp([name,'.Rxx =']);
fprintf("\n");
Rmxx = max([size(Rxxstr_o,2),size(Rxxstr_a,2),size(Rxxstr_b,2)]);
if size(Rxxstr_o,2)<Rmxx
    Rxxstr_o(end+1:Rmxx) = " ";
end
if size(Rxxstr_a,2)<Rmxx
    Rxxstr_a(end+1:Rmxx) = " ";
end
if size(Rxxstr_b,2)<Rmxx
    Rxxstr_b(end+1:Rmxx) = " ";
end
Rxxrowmax_o = max(arrayfun(@(x) strlength(x),Rxxstr_o));
Rxxrowmax_a = max(arrayfun(@(x) strlength(x),Rxxstr_a));
Rxxrowmax_b = max(arrayfun(@(x) strlength(x),Rxxstr_b));
for i = 1:Rmxx
    fprintf("%3s %*s | %*s | %*s \n", " ", Rxxrowmax_o, Rxxstr_o(i), Rxxrowmax_a, Rxxstr_a(i), Rxxrowmax_b, Rxxstr_b(i));
end

%%% Rx2
fprintf("\n");
disp([name,'.Rx2 =']);
fprintf("\n");
Rmx2 = max([size(Rx2str_o,2),size(Rx2str_a,2),size(Rx2str_b,2)]);
if size(Rx2str_o,2)<Rmx2
    Rx2str_o(end+1:Rmx2) = " ";
end
if size(Rx2str_a,2)<Rmx2
    Rx2str_a(end+1:Rmx2) = " ";
end
if size(Rx2str_b,2)<Rmx2
    Rx2str_b(end+1:Rmx2) = " ";
end
Rx2rowmax_o = max(arrayfun(@(x) strlength(x),Rx2str_o));
Rx2rowmax_a = max(arrayfun(@(x) strlength(x),Rx2str_a));
Rx2rowmax_b = max(arrayfun(@(x) strlength(x),Rx2str_b));
for i = 1:Rmx2
    fprintf("%3s %*s | %*s | %*s \n", " ", Rx2rowmax_o, Rx2str_o(i), Rx2rowmax_a, Rx2str_a(i), Rx2rowmax_b, Rx2str_b(i));
end

%%% Ryy
fprintf("\n");
disp([name,'.Ryy =']);
fprintf("\n");
Rmyy = max([size(Ryystr_o,2),size(Ryystr_a,2),size(Ryystr_b,2)]);
if size(Ryystr_o,2)<Rmyy
    Ryystr_o(end+1:Rmyy) = " ";
end
if size(Ryystr_a,2)<Rmyy
    Ryystr_a(end+1:Rmyy) = " ";
end
if size(Ryystr_b,2)<Rmyy
    Ryystr_b(end+1:Rmyy) = " ";
end
Ryyrowmax_o = max(arrayfun(@(x) strlength(x),Ryystr_o));
Ryyrowmax_a = max(arrayfun(@(x) strlength(x),Ryystr_a));
Ryyrowmax_b = max(arrayfun(@(x) strlength(x),Ryystr_b));
for i = 1:Rmyy
    fprintf("%3s %*s | %*s | %*s \n", " ", Ryyrowmax_o, Ryystr_o(i), Ryyrowmax_a, Ryystr_a(i), Ryyrowmax_b, Ryystr_b(i));
end

%%% Ry2
fprintf("\n");
disp([name,'.Ry2 =']);
fprintf("\n");
Rmy2 = max([size(Ry2str_o,2),size(Ry2str_a,2),size(Ry2str_b,2)]);
if size(Ry2str_o,2)<Rmy2
    Ry2str_o(end+1:Rmy2) = " ";
end
if size(Ry2str_a,2)<Rmy2
    Ry2str_a(end+1:Rmy2) = " ";
end
if size(Ry2str_b,2)<Rmy2
    Ry2str_b(end+1:Rmy2) = " ";
end
Ry2rowmax_o = max(arrayfun(@(x) strlength(x),Ry2str_o));
Ry2rowmax_a = max(arrayfun(@(x) strlength(x),Ry2str_a));
Ry2rowmax_b = max(arrayfun(@(x) strlength(x),Ry2str_b));
for i = 1:Rmy2
    fprintf("%3s %*s | %*s | %*s \n", " ", Ry2rowmax_o, Ry2str_o(i), Ry2rowmax_a, Ry2str_a(i), Ry2rowmax_b, Ry2str_b(i));
end

%%% R2x
fprintf("\n");
disp([name,'.R2x =']);
fprintf("\n");
Rm2x = max([size(R2xstr_o,2),size(R2xstr_a,2),size(R2xstr_b,2)]);
if size(R2xstr_o,2)<Rm2x
    R2xstr_o(end+1:Rm2x) = " ";
end
if size(R2xstr_a,2)<Rm2x
    R2xstr_a(end+1:Rm2x) = " ";
end
if size(R2xstr_b,2)<Rm2x
    R2xstr_b(end+1:Rm2x) = " ";
end
R2xrowmax_o = max(arrayfun(@(x) strlength(x),R2xstr_o));
R2xrowmax_a = max(arrayfun(@(x) strlength(x),R2xstr_a));
R2xrowmax_b = max(arrayfun(@(x) strlength(x),R2xstr_b));
for i = 1:Rm2x
    fprintf("%3s %*s | %*s | %*s \n", " ", R2xrowmax_o, R2xstr_o(i), R2xrowmax_a, R2xstr_a(i), R2xrowmax_b, R2xstr_b(i));
end

%%% R2y
fprintf("\n");
disp([name,'.R2y =']);
fprintf("\n");
Rm2y = max([size(R2ystr_o,2),size(R2ystr_a,2),size(R2ystr_b,2)]);
if size(R2ystr_o,2)<Rm2y
    R2ystr_o(end+1:Rm2y) = " ";
end
if size(R2ystr_a,2)<Rm2y
    R2ystr_a(end+1:Rm2y) = " ";
end
if size(R2ystr_b,2)<Rm2y
    R2ystr_b(end+1:Rm2y) = " ";
end
R2yrowmax_o = max(arrayfun(@(x) strlength(x),R2ystr_o));
R2yrowmax_a = max(arrayfun(@(x) strlength(x),R2ystr_a));
R2yrowmax_b = max(arrayfun(@(x) strlength(x),R2ystr_b));
for i = 1:Rm2y
    fprintf("%3s %*s | %*s | %*s \n", " ", R2yrowmax_o, R2ystr_o(i), R2yrowmax_a, R2ystr_a(i), R2yrowmax_b, R2ystr_b(i));
end

%%% R22
fprintf("\n");
disp([name,'.R22 =']);
fprintf("\n");
Rm22 = 0;
R22rowmax = zeros(3,3);
for i=1:3
    for j=1:3
        Rm22 = max(Rm22,size(R22str{i,j},2));
        R22rowmax(i,j) = max(arrayfun(@(x) strlength(x),R22str{i,j}));
    end
end
for i=1:3
    for j=1:3
        if size(R22str{i,j},2)<Rm22
            Rstr_ij = R22str{i,j};
            Rstr_ij(end:1:Rm22) = " ";
            R22str{i,j} = Rstr_ij;
        end
    end
end
R22rowmax = max(R22rowmax,[],1);

divlen = sum(R22rowmax);
div = "";
for i=1:divlen+2*3
div = div+"-";
end

% Display the first row
for k=1:Rm22
    R22strj = R22str{1,1};
    fprintf("%3s %*s ", " ", R22rowmax(1), R22strj(k));
    for j=2:3
        R22strj = R22str{1,j};
        fprintf("| %*s ", R22rowmax(j), R22strj(k));
    end
    fprintf("\n");
end
% Display the remaining rows
for i=2:3
    fprintf("%3s %s \n"," ",div); 
    for k = 1:Rm22
        R22strj = R22str{i,1};
        fprintf("%3s %*s ", " ", R22rowmax(1), R22strj(k));
        for j=2:3
            R22strj = R22str{i,j};
            fprintf("| %*s ", R22rowmax(j), R22strj(k));
        end
        fprintf("\n");
    end
end


end



function [str,big] = getstringout(prop,big)
llimit = 50; %80
if isempty(prop)
str = "[]"; 
else
if ~big
    for i=1:size(prop,1)
        row ="";
        for j=1:size(prop,2)
            temp = evalc('prop(i,j)');
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
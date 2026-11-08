function show(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show(P,opts) function displays components specified in opts
% of the opvar object P
% 
% INPUT
% P : An opvar variable
% -Optional INPUTS
% opts: flag to decide which components to be printed (default selection -
% everything)
% 
% OUTPUT
% command line display of elements of opvar
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - show
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
% Initial coding MMP, SS  - 6_26_2020

if nargin==1
opts = {'P','Q1','Q2','R0','R1','R2'};
end

if ~isa(P,'opvar')
    error("Input object must be of opvar class");
end

name = inputname(1);
for i=opts
    if contains(i,'R')
        disp([name,'.R.',i{:},'=']);
        P.R.(i{:})
    else
        disp([name,'.',i{:},'=']);
        P.(i{:})
    end
end
end
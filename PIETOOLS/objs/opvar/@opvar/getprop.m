function out=getprop(T,prop)
% out=get(T,prop) retrieves the components of the opvar object T
% 
% INPUTS:
% T : An opvar variable
% prop: structure with property to be accessed
% 
% OUTPUTS:
% out: property value
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

if length(prop)==1
    prop = prop.subs;
    indx = 1:2; indy = 1:2;
elseif length(prop)==2
    if strcmp(prop(1).subs,'dim')
        indx = prop(2).subs{1}; indy = prop(2).subs{2};
        prop = 'dim';
    elseif prop(1).subs=='I' && strcmp(prop(2).type,'()')
        prop = [prop(1).subs,'(',num2str(prop(2).subs{1}),')'];
    else
        prop = [prop(1).subs '.' prop(2).subs];
    end
else
    error('Invalid property name');
end

switch prop
    case 'P'
        out = T.P;
    case 'Q1'
        out = T.Q1;
    case 'Q2'
        out = T.Q2;
    case 'R'
        out = T.R;
    case 'R.R0'
        out = T.R.R0;
    case 'R.R1'
        out = T.R.R1;
    case 'R.R2'
        out = T.R.R2;
    case 'dim'
        out = T.dim(indx,indy);
    case 'I'
        out = T.I;
    case 'I(1)'
        out = T.I(1);
    case 'I(2)'
        out = T.I(2);
    case 'var1'
        out = T.var1;
    case 'var2'
        out = T.var2;
    case 'dimdependent'
        out = T.dimdependent;
    otherwise
        error('Invalid property name');
end
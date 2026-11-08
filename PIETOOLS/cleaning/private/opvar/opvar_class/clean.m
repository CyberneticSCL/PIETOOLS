function out = clean(in,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = clean(in,tol) removes all polynomial terms with a coefficient smaller than tol
% 
% INPUT
% in: opvar class object
% tol: tolerance value, default 1e-7
% 
% OUTPUT
% out: opvar class object
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PIETools - isvalid
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

if nargin<2
    tol=1e-7;
end
out = in;
fields = {'P','Q1','Q2','R0','R1','R2'};

for i=fields
    if any(ismember(i{:},'R'))
        tmp = in.R.(i{:});
    else
        tmp = in.(i{:});
    end
    if isa(tmp,'double')
        tmp(find(abs(tmp)<tol))=0;
    else
        tmp.coefficient(find(abs(tmp.coefficient)<tol)) = 0;
    end
    if any(ismember(i{:},'R'))
        out.R.(i{:})=tmp;
    else
        out.(i{:})=tmp;
    end
end
end
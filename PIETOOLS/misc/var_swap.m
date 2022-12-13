function b = var_swap(a,old,new)
% function B = subs(A,Old,New);
%
% DESCRIPTION 
%  Swap two symbolic variables in polynomial expression.
%   
% INPUTS 
%   A: polynomial
%   Old: single polynomial variable
%   New: single polynomial variable
%
% OUTPUTS  
%   B:  polynomial
%       B(old,new)=A(new,old)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% 6/23/2013: MMP Coding.
% 09/26/2021: Added call to dpvar version

  if isempty(a)
      b=a;
      return
  end
  
  if isa(a,'dpvar')
      b = varswap(a,old,new);
      return
  end
      
  a = polynomial(a);

  temp_varname=a.varname;
  old_idx=find(strcmp(old.varname,temp_varname));  
  new_idx=find(strcmp(new.varname,temp_varname));  
  temp_varname(new_idx)=cellstr(old.varname);
  temp_varname(old_idx)=cellstr(new.varname);
  %temp(lind1)=cellstr([repmat('coeff_',length(lind1),1), int2str(sos.var.idx{var}-sos.var.idx{1}+lind1)]);
  b=polynomial(a.coefficient,a.degmat,temp_varname,a.matdim);
  

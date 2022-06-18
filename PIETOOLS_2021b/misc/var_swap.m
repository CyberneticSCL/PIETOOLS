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
  

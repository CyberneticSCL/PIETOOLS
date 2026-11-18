function n = numArgumentsFromSubscript(obj,s,indexingContext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This method is automatically called when one of the following statements
% are executed:
% a) obj.s — Number of elements referenced in a statement (indexingContext = Statement) 
% b) (obj.s) — Number of elements returned in an expression (indexingContext = Expression) 
% c) [obj.s] = rhs — Number of values assigned with a comma-separated list (indexingContext = Assignment)

   n=1;
end 
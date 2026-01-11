function Cop = minus(Aop,Bop)
% COP = MINUS(AOP,BOP) returns the ndopvar object representing the
% difference of the PI operators defined by AOP and BOP.

Cop = plus(Aop,uminus(Bop));

end


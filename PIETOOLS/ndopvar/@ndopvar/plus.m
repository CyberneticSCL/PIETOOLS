function Cop = plus(Aop,Bop)
% COP = PLUS(AOP,BOP) returns the 'ndopvar' object COP representing the sum
% of the PI operators defined by 'ndopvar' objects AOP and BOP


% Check that the operators can indeed be added
if any(Aop.dim~=Bop.dim)
    error("Dimensions of the operators must match.")
end
if size(Aop.dom,1)~=size(Bop.dom,1) || any(any(Aop.dom~=Bop.dom))
    error("Spatial domains on which operators are defined should match.")
end

if any(Aop.deg~=Bop.deg)
    error("Addition of operators with different monomial degrees is currently not supported.")
elseif ~isempty(Aop.dvarname) || ~isempty(Bop.dvarname)
    error("Addition of operators with decision variables is currently not supported.")
end

% Declare the sum
Cop = Aop;
for ii=1:numel(Cop.C)
    Cop.C{ii} = Aop.C{ii} + Bop.C{ii};
end

end
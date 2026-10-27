%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_3PI_opmult_discretize.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs A - discrete matrix representations of a polymonial
% multiplicative operator component of 3PI operator
% Discretizes an operator component (NOT FULL OPERATOR)
%
% Inputs:
% N   - polynomial order of Chebyshev discretization polynomial
% Rop -  polynomial multplicative operator R0 of size 1x1
%
% Outputs:
% A - square matrix of size N+1 x N+1 that represents a Chebyshev
% discretizaiton of the polynomial multiplicative matrix operator. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 06_28_2022
% YP -  added support for an arbitary variable name - 04_16_2024
% DJ, 12/16/2024: Make sure variable to substitute matches that of Rop;
% YP, 2/1/2026: changed to square operator, truncation done on assembled
% operators

function A=PIESIM_3PI_opmult_discretize(N, Rop)

A(1:N+1,1:N+1)=0;

if isa(Rop,'polynomial')
    deg=Rop.maxdeg;
else
    deg=1;
end

chebgrid=cos(pi*(0:deg+1)/(deg+1));
if isa(Rop,'polynomial') && ~isdouble(Rop)
    if numel(Rop.varname)>1                                                 % DJ, 12/16/2024
        error("The input polynomial depends on multiple variables; discretization not supported...")
    else
        var = polynomial(Rop.varname);
        Reval = subs(Rop,var,chebgrid);
    end    
else
    Reval=Rop*ones(1,deg+2);
end

acheb=fcht(double(Reval));

for i=1:N+1
    id=i-1;
    for j=1:length(acheb)
        jd=j-1;
        if (jd+id<=N)
        A(jd+id+1,i)=A(jd+id+1,i)+0.5*acheb(j);
        end
        if (abs(jd-id)<=N)
        A(abs(jd-id)+1,i)=A(abs(jd-id)+1,i)+0.5*acheb(j);
        end
    end
end





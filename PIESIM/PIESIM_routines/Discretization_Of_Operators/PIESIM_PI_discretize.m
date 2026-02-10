%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_PI_discretize.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs Rcheb - PI operator evaluated on a Chebyshev grid
%
% Inputs:
% R - compnent PI operator to be discretized
% maxdeg - maximum polynomial degree of a PI operator R
%
% Output:
% Rcheb - a square matrix of size (maxdeg+2)x(maxdeg+2) that evaluates the PI operator
% on a Chebysheb grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Used by PIESIM_3PI2Mat_cheb_opint_discretize.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 12_21_2024
function Rcheb=PIESIM_PI_discretize(R,maxdeg)

chebgrid=cos(pi*(0:maxdeg+1)/(maxdeg+1));
if isa(R,'polynomial')
    switch(R.nvars)
        case 2
        % Two variables present
            Reval=subs(subs(R,R.varname(2),chebgrid),R.varname(1),chebgrid);
        case 1
        % One variable present
            Re=subs(R,R.varname,chebgrid);
            var=cell2mat(R.varname);
            if (length(var)<=2)
                Reval=repmat(Re',1,maxdeg+2);
            else
                Reval=repmat(Re,maxdeg+2,1);
            end
        case 0 
            % Scalar polynomial
            Reval=R*ones(maxdeg+2);
end
else
    Reval=R*ones(maxdeg+2);
end

Reval=double(Reval);
Rcheb=fcgltran2d(Reval,1);
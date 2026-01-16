function Pop = scalar(alpha,Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pop = scalar(alpha,Pop) returns the 'nopvar' object Pop representing 
% the scalar product alpha*Pop for the PI operator defined by Pop.
% Version: 1.0
% 
% INPUT
% alpha: double to scale coefficients in Pop.C;
% Pop:   nopvar object;
% 
% OUTPUT
% Pop:   scaled nopvar object;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding CR - 1/15/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:numel(Pop.C)
    Pop.C{ii} = alpha*Pop.C{ii};
end

end
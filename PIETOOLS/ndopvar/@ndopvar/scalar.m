function Pop = scalar(alpha,Pop) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pop = scalar(alpha,Pop) returns the 'nopvar' object representing 
% the scalar product alpha*Pop for the PI operator defined by Pop.
% Date: 1/15/26
% Version: 1.0
% 
% INPUT
% alpha: double
% Pop:   nopvar object
% 
% OUTPUT
% Pop:   nopvar object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pop.C = alpha*Pop.C;

end
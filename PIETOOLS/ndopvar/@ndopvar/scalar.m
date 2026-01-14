function Pop = scalar(alpha,Pop) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pop = scalar(alpha,Pop) returns the 'ndopvar' object representing 
% the scalar product alpha*Pop for the PI operator defined by Pop.
% Date: 1/14/26
% Version: 1.0
% 
% INPUT
% alpha: double
% Pop:   ndopvar object
% 
% OUTPUT
% Pop:   ndopvar object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:numel(Pop.C)
    Pop.C{ii} = alpha*Pop.C{ii};
end
end
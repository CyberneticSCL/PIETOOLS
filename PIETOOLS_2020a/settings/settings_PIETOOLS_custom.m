%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_order1 = 1; 
n_order2 = 1; % This is supposed to be an accuracy/non-balancing degree
n_order3 = 1;
n_order4 = max([n_order2,n_order3]); n_order4=4;


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for the LF

override1=1;   % No Psatz in LF
options1.sep = 0; %this is to select separable case, R1=R2
options1.exclude=[0 0 0 0]; % using LU or UL decomposition in the LF seems to have minimal impact...
options12.sep = 0; %this is to select separable case, R1=R2
options12.exclude = [0 0 0 0]; % items to exclude from the PSatz part of LF
options12.psatz=1; % this should always be 1, otherwise, what is the point?
dd1 = {n_order1, [n_order2, n_order3, n_order4],[n_order2, n_order3, n_order4]}; %use above array format if needed- sachin
%dd12 = {n_order1-1, [n_order2-1, n_order3-1, n_order4-1],[n_order2-1, n_order3-1, n_order4-1]}; %use above array format if needed- sachin
dd12 = {n_order1, [n_order2, n_order3, n_order4],[n_order2, n_order3, n_order4]}; %use above array format if needed- sachin

 

% new files need exclude instead of full and pure - sach

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for the indeterminate variable, Z
% This is the degree choices for the indefinite variable Z used in
% controller and estimator synthesis
ddZ=[2*n_order1 2*n_order2 2*n_order3];
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for the derivative
%
sosineq_on=0;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sosineq_on
    % These options passed to sosineq
    opts.psatz = 0;
    opts.pure = 1;
    % if sosineq is used, then options2, options3 defined here are unused
else
    % These options are NOT passed to sosineq
    Dup=5;
    override2=0;   % =0 to include Psatz in Derv
    % when not using sosineq, options 2 is for the nominal part of the derv,
    % options3 is for the Psatz part of the derv
    options2.psatz = 0;
    options2.exclude= [0 1 0 0];
    options2.sep= 0;
    dd2 = {n_order1+Dup, [n_order2+Dup, n_order3+Dup, n_order4+Dup], [n_order2+Dup, n_order3+Dup, n_order4+Dup]};
    
    options3.psatz=1; % this should always be 1, otherwise no point
    options3.exclude = [0 1 0 0]; % using the LU or UL decomposition in the psatz term seems to have minimal impact when norder2=1 norder1=2
   dd3 = {n_order1+Dup, [n_order2+Dup, n_order3+Dup, n_order4+Dup], [n_order2+Dup, n_order3+Dup, n_order4+Dup]};
%     dd3 = {n_order1+Dup-1, [n_order2+Dup-1, n_order3+Dup-1, n_order4+Dup-1], [n_order2+Dup-1, n_order3+Dup-1, n_order4+Dup-1]};
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

dd1 = {1}; dd12 = {1}; dd2 = {1,5}; dd3 = {1,5};



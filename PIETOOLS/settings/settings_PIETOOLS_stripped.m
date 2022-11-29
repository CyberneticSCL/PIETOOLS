function settings = settings_PIETOOLS_stripped()
% Returns settings structure with degrees between 2 and 3 with lower
% triangular structure in square root of positive definite PIs
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
% DJ - 06/02/2021 - adjusted to store sosineq_on and clear auxilliary variables at the end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_order1 = 1; 
n_order2 = 1; % This is supposed to be an accuracy/non-balancing degree
n_order3 = 0;
n_order4 = max([n_order2,n_order3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for the LF

settings.override1=1;   % No Psatz in LF
settings.options1.sep = 0; %this is to select separable case, R1=R2
settings.options1.exclude=[0 0 0 0]; % using LU or UL decomposition in the LF 
settings.options12.sep = 0; %this is to select separable case, R1=R2
settings.options12.exclude = [0 0 0 0]; % items to exclude from the PSatz part of LF
settings.options12.psatz=1; % this should always be 1, otherwise, what is the point?
settings.dd1 = {n_order1, [n_order2, n_order3, n_order4],[n_order2, n_order3, n_order4]}; 
% dd12 = {n_order1-1, [n_order2-1, n_order3-1, n_order4-1],[n_order2-1, n_order3-1, n_order4-1]}; 
settings.dd12 = {n_order1, [n_order2, n_order3, n_order4],[n_order2, n_order3, n_order4]}; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for the indeterminate variable, Z
% This is the degree choices for the indefinite variable Z used in
% controller and estimator synthesis
settings.ddZ=[2*n_order1 2*n_order2 2*n_order3];
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for the derivative
%
if ~exist('sosineq_on','var')       
    sosineq_on = 0;
end
settings.sosineq_on = sosineq_on;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sosineq_on
    % These options passed to sosineq
    settings.opts.psatz = 0;
    settings.opts.pure = 1;
    % if sosineq is used, then options2, options3 defined here are unused
else
    % These options are NOT passed to sosineq
    Dup=1;
    settings.override2=1;   % Psatz in Derv
    % when not using sosineq, options 2 is for the nominal part of the derv,
    % options3 is for the Psatz part of the derv
    settings.options2.exclude= [0 1 0 0];
    settings.options2.sep= 0;
    settings.dd2 = {n_order1+Dup, [n_order2+Dup-1, n_order3+Dup, n_order4+Dup], [n_order2+Dup-1, n_order3+Dup, n_order4+Dup]};
    
    settings.options3.psatz=1; % this should always be 1, otherwise no point
    settings.options3.exclude = [0 1 0 1]; % using the LU or UL decomposition in the psatz term seems to have minimal impact when norder2=1 norder1=2
    settings.dd3 = {n_order1+Dup-1, [n_order2+Dup-2, n_order3+Dup-1, n_order4+Dup-1], [n_order2+Dup-2, n_order3+Dup-1, n_order4+Dup-1]};
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Clear internal variables
clear Dup n_order1 n_order2 n_order3 n_order4



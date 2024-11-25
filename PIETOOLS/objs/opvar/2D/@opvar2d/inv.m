function Pinv = inv(Pop,inv_tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pinv] = inv(Pop) takes in one operator, 
% Pop: R^m0 x L2^mx x L2^my x L2^m2   to    R^n0 x L2^nx x L2^ny x L2^n2
% It returns the (left) inverse 
% Pinv: R^n0 x L2^nx x L2^ny x L2^n2   to    R^m0 x L2^mx x L2^my x L2^m2
% 
% INPUT
% P: opvar2d class object
%
% OUTPUT
% Pinv: the matlab structure such that subsequent application Pinv*P*y for
% an arbitrary polynomial y:R^n0 x L2^nx(s1) x L2^ny(s2) returns y.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - inv
%
% Copyright (C)2021  M. Peet, D. Jagt
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
%
% Initial coding DJ - 07_01_2024


% % Deal with empty operator...
if isempty(Pop)
    Pinv = Pop;
    return
end
% % Set a tolerance if not specified.
if nargin==1
    inv_tol = 1e-8;
end

% % % Check if the operator is separable.
% % % If so, we can compute an exact inverse (if the diagonal terms are
% % % invertible...)
if is_separable(Pop)
    % % If the operator is separable, use the appropriate function.
    Pinv = inv_opvar2d_separable(Pop,inv_tol);
else
    % % Otherwise, see if we can get an inverse using mrdivide.
    deg_fctr = 2*ones(4,2);
    deg_fctr_max = [10,10; 10,10; 5,5; 5,5];
    Pinv = mrdivide(1,P1op,deg_fctr,inv_tol,deg_fctr_max);
end

end
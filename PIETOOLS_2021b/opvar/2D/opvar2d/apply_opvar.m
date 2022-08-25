function Px = apply_opvar(P,v)
% Px = apply_opvar(P,x) applies the PI operator P to polynomial x and
% returns the result Px=P*x;
%
% INPUT
% - P:      m x n opvar or opvar2d class object.
% - v:      n x k 'polynomial' or 'double' class object.
%
% OUTPUT:
% - Pv:     m x k 'polynomial' class object, corresponding to the result of
%           applying the PI operator associated to P to the polynomial x.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - apply_opvar
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/23/2022
%

% Check that the PI operator is properly specified, and whether it is a 1D
% or 2D operator.
use_2d = false;
if ~isa(P,'opvar') && ~isa(P,'opvar2d')
    error('The first argument should be of class ''opvar'' or ''opvar2d'', representing a PI operator.')
elseif isa(P,'opvar2d')
    use_2d = true;
end

% Check that the polynomial has been properly specified.
if isa(v,'double')
    v = polynomial(v);
elseif ~isa(v,'polynomial')
    error('The second argument should be of class ''double'' or ''polynomial'', representing a polynomial on which to apply the PI operator.')
end

% Check that the dimensions of the operator and polynomial match.
nc_arr = P.dim(:,2);
nnc_arr = cumsum(nc_arr);
if size(v,1)~=nnc_arr(end)
    error('The column dimension of the input PI operator should match the number of rows of the polynomial.')
end
if all(nnc_arr==0)
    nr_arr = P.dim(:,1);
    Px = polynomial(zeros(sum(nr_arr),size(v,2)));
    return
end

% Distinguish 1D and 2D cases.
if ~use_2d
    % Split the polynomial into finite-dimensional and 1D components.
    v0 = v(1:nnc_arr(1),1);
    vx = v(nnc_arr(1)+1:nnc_arr(2),1); 
    
    % Extract the spatial variables and their domain.
    ds1 = P.var1;   dt1 = P.var2;
    a1 = P.I(1);    b1 = P.I(2);
    
    % % Apply the operator:
    % Finite-dimensional output.
    Px0 = P.R00*v0 + int(P.R0x*vx,ds1,a1,b1);

    % Infinite-dimensional (1D) output.
    Pxx = P.Rx0*v0;
    Pxx = Pxx + P.Rxx{1,1}*vx + int(P.Rxx{2,1}*subs(vx,ds1,dt1),dt1,a1,ds1) + int(P.Rxx{3,1}*subs(vx,ds1,dt1),dt1,ds1,b1);
    
    % Return the result.
    Px = [Px0;Pxx];
    
else
    % Split the polynomial into finite-dimensional, 1D, and 2D components
    v0 = v(1:nnc_arr(1),1);
    vx = v(nnc_arr(1)+1:nnc_arr(2),1);
    vy = v(nnc_arr(2)+1:nnc_arr(3),1);
    v2 = v(nnc_arr(3)+1:nnc_arr(4),1);

    % Extract the spatial variables
    ds = P.var1;        dt = P.var2;
    ds1 = ds(1);        ds2 = ds(2);
    dt1 = dt(1);        dt2 = dt(2);

    % Establish the lower and upper boundaries of the domain of the variables.
    a1 = P.I(1,1);      a2 = P.I(2,1);
    b1 = P.I(1,2);      b2 = P.I(2,2);

    % % Apply the operator:
    % Finite-dimensional output.
    Px0 = P.R00*v0 + int(P.R0x*vx,ds1,a1,b1) + int(P.R0y*vy,ds2,a2,b2) + int(int(P.R02*v2,ds1,a1,b1),ds2,a2,b2);

    % 1D output in first variable.
    Pxx = P.Rx0*v0 + int(P.Rxy*vy,ds2,a2,b2);
    Pxx = Pxx + P.Rxx{1,1}*vx + int(P.Rxx{2,1}*subs(vx,ds1,dt1),dt1,a1,ds1) + int(P.Rxx{3,1}*subs(vx,ds1,dt1),dt1,ds1,b1);
    Pxx = Pxx + int(P.Rx2{1,1}*v2 + int(P.Rx2{2,1}*subs(v2,ds1,dt1),dt1,a1,ds1) + int(P.Rx2{3,1}*subs(v2,ds1,dt1),dt1,ds1,b1),ds2,a2,b2);

    % 1D output in second variable.
    Pxy = P.Ry0*v0 + int(P.Ryx*vx,ds1,a1,b1);
    Pxy = Pxy + P.Ryy{1,1}*vy + int(P.Ryy{1,2}*subs(vy,ds2,dt2),dt2,a2,ds2) + int(P.Ryy{1,3}*subs(vy,ds2,dt2),dt2,ds2,b2);
    Pxy = Pxy + int(P.Ry2{1,1}*v2 + int(P.Ry2{1,2}*subs(v2,ds2,dt2),dt2,a2,ds2) + int(P.Ry2{1,3}*subs(v2,ds2,dt2),dt2,ds2,b2),ds1,a1,b1);
    
    % 2D output.
    Px2 = P.R20*v0;
    Px2 = Px2 + P.R2x{1,1}*vx + int(P.R2x{2,1}*subs(vx,ds1,dt1),dt1,a1,ds1) + int(P.R2x{3,1}*subs(vx,ds1,dt1),dt1,ds1,b1);
    Px2 = Px2 + P.R2y{1,1}*vy + int(P.R2y{1,2}*subs(vy,ds2,dt2),dt2,a2,ds2) + int(P.R2y{1,3}*subs(vy,ds2,dt2),dt2,ds2,b2);
    Px2 = Px2 + P.R22{1,1}*v2 + int(P.R22{2,1}*subs(v2,ds1,dt1),dt1,a1,ds1) + int(P.R22{3,1}*subs(v2,ds1,dt1),dt1,ds1,b1) ...
        + int(P.R22{1,2}*subs(v2,ds2,dt2),dt2,a2,ds2) + int(P.R22{1,3}*subs(v2,ds2,dt2),dt2,ds2,b2) ...
        + int(int(P.R22{2,2}*subs(v2,ds,dt),dt1,a1,ds1),dt2,a2,ds2) + int(int(P.R22{3,2}*subs(v2,ds,dt),dt1,ds1,b1),dt2,a2,ds2) ...
        + int(int(P.R22{2,3}*subs(v2,ds,dt),dt1,a1,ds1),dt2,ds2,b2) + int(int(P.R22{3,3}*subs(v2,ds,dt),dt1,ds1,b1),dt2,ds2,b2);

    % Return the result.
    Px = [Px0;Pxx;Pxy;Px2];
end

end
function Pop_out = opvar2d2multiplier(Pop_in,tol)
% OPVAR2D2MULTIPLIER converts an opvar2d object representing a PI operator
%   Pop_in: RxL2[x]xL2[y]xL2[x,y] --> RxL2[x]xL2[y]xL2[x,y] 
% that acts only as a multiplier to an associated PI operator
%   Pop_out:  R^4 --> RxL2[x]xL2[y]xL2[x,y]
% that still acts only as a multiplier, but now mapping finite-dimensional
% vectors to functions.
%   
% INPUT
% - Pop_in:     opvar2d class object representing a PI operator
%                   RxL2[x]xL2[y]xL2[x,y] --> RxL2[x]xL2[y]xL2[x,y] 
%               that acts only as a multiplier operator. That is, only
%               parameters Pop_in.R00, Pop_in.Rx0, Pop_in.Rxx{1},
%               Pop_in.Ry0, Pop_in.Ryy{1}, Pop_in.R20, Pop_in.R2x{1},
%               Pop_in.R2y{1} and Pop_in.R22{1,1} may be nonzero, as all
%               other parameters correspond to integral operators.
% - tol:        Scalar double specifying the tolerance below which
%               parameters are assumed to be zero. Defaults to 1e-12.
%
% OUTPUT:
% - Pop_out:    opvar2d class object representing a PI operator
%                   R --> RxL2[x]xL2[y]xL2[x,y] 
%               where all the multipliers in Pop_in have been concatenated
%               column-wise into just the fields
%                   Pop_out.R00 = Pop_in.R00;
%                   Pop_out.Rx0 = [Pop_in.Rx0, Pop_in.Rxx{1}];
%                   Pop_out.Ry0 = [Pop_in.Ry0, Pop_in.Ryy{1}];
%                   Pop_out.R20 = [Pop_in.R20, Pop_in.R2x{1},...
%                                       Pop_in.R2y{1}, Pop_in.R22{1,1}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - opvar2d2multiplier
%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 09/07/2024
%

% % Process the inputs.
if ~isa(Pop_in,'opvar2d')
    error("Input must be of type 'opvar2d'.")
end
if nargin<=1
    % Set default tolerance.
    tol = 1e-12;
end
if ~any(Pop_in.dim(:,1)) || ~any(Pop_in.dim(2:end,2))
    % Deal with operator that is empty or already maps from \R.
    Pop_out = Pop_in;
    return
end

% % Check that the operator acts only as a multiplier up to tolerance.
% Get rid of everything below tolerance.
Pop_in = clean_opvar(Pop_in,tol);

% Store all parts that act as multiplier in a temporary operator
Pop_tmp = opvar2d();
Pop_tmp.I = Pop_in.I;       Pop_tmp.dim = Pop_in.dim;
Pop_tmp.var1 = Pop_in.var1; Pop_tmp.var2 = Pop_in.var2;

Pop_tmp.R00 = Pop_in.R00;
Pop_tmp.Rx0 = Pop_in.Rx0;           Pop_tmp.Rxx{1} = Pop_in.Rxx{1};
Pop_tmp.Ry0 = Pop_in.Ry0;           Pop_tmp.Ryy{1} = Pop_in.Ryy{1};
Pop_tmp.R20 = Pop_in.R20;           Pop_tmp.R22{1,1} = Pop_in.R22{1,1};
Pop_tmp.R2x{1} = Pop_in.R2x{1};     Pop_tmp.R2y{1} = Pop_in.R2y{1};

% Check that the multiplier components indeed comprise the full operator.
if ~(Pop_tmp==Pop_in)
    error('Specified operator acts as integral; conversion to multiplier operator is not supported.')
end


% % Build the multiplier operator associated to Pop_in
% Initialize an operator mapping \R --> \R x L2[x] x L2[y] x L2[x,y]
nr_op = Pop_in.dim(:,1);        nc_op = Pop_in.dim(:,2);
Pop_out = opvar2d();
Pop_out.I = Pop_in.I;       
Pop_out.dim = [nr_op,[sum(nc_op);0;0;0]];
Pop_out.var1 = Pop_in.var1; 
Pop_out.var2 = Pop_in.var2;

% Move all parameters acting as multiplier to the left columns
Pop_out.R00 = [Pop_in.R00, zeros(nr_op(1),sum(nc_op(2:end)))];
Pop_out.Rx0 = [Pop_in.Rx0, Pop_in.Rxx{1}, zeros(nr_op(2),sum(nc_op(3:end)))];
Pop_out.Ry0 = [Pop_in.Ry0, zeros(nr_op(3),nc_op(2)), Pop_in.Ryy{1}, zeros(nr_op(3),nc_op(4))];
Pop_out.R20 = [Pop_in.R20, Pop_in.R2x{1}, Pop_in.R2y{1}, Pop_in.R22{1,1}];

end
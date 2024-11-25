function [is_sep,is_full_int,is_no_mult] = is_separable(Pop,Ctol)
% IS_SEP = IS_SEPARABLE(POP,CTOL) Checks if an 'opvar2d' object "Pop"
% defines a separable 2D PI operator.
%
% INPUT
% - Pop:    'opvar2d' class object, representing a 2D PI operator.
% - Ctol:   A positive scalar value, indicating a "tolerance" for nonzero
%           coefficients. For any parameter in 'Pop', if the coefficients
%           defining this parameter are smaller than 'Ctol', then they are
%           interpreted to be zero. Defaults to 1e-12.
%
% OUTPUT
% - is_sep: Boolean variable specifying whether the input object 'Pop'
%           represents a separable 2D PI operator. This is the case if each
%           sub-component RRij: L2[Dom_in] --> L2[Dom_out] takes the form
%             (Rij*u)(x) = Fij(x)*int_{Dom_in}Gij(y)*u(y)dy    x in Dom_out
%           Diagonal components may involve a multiplier operator, as
%             (Rii*u)(x) = Rii0(x)*u(x) +Fij(x)*int_{Dom_in}Gij(x)*u(x)dx.
%
% NOTES
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
%
% Initial coding DJ - 07/01/2024


% % Convert all operator parameters to class 'polynomial'
Pop = poly_opvar2d(Pop);

% Coefficients smaller than Ctol are assumed to be zero.
if nargin==1
    Ctol = 1e-12;
end

% % Check that the operator only involves full integrals.
is_full_int = true(5,1);
if ~all(all(isequal(cleanpoly(Pop.Rxx{2}-Pop.Rxx{3},Ctol),0)))
    % Check if the map L2[x]-->L2[x] is defined by full integral operator
    is_full_int(1) = false;
elseif ~all(all(isequal(cleanpoly(Pop.Ryy{2}-Pop.Ryy{3},Ctol),0)))
    % Check if the map L2[y]-->L2[y] is defined by full integral operator
    is_full_int(2) = false;
elseif ~all(all(isequal(cleanpoly(Pop.R22{2,2}-Pop.R22{3,2},Ctol),0))) ||...
        ~all(all(isequal(cleanpoly(Pop.R22{2,2}-Pop.R22{2,3},Ctol),0))) ||...
            ~all(all(isequal(cleanpoly(Pop.R22{2,2}-Pop.R22{3,3},Ctol),0)))
    % Check if the map L2[x,y]-->L2[x,y] is defined by full integral operator
    is_full_int(3) = false;
elseif ~all(all(isequal(cleanpoly(Pop.Rx2{2}-Pop.Rx2{3},Ctol),0))) ||...
        ~all(all(isequal(cleanpoly(Pop.R2x{2}-Pop.R2x{3},Ctol),0)))
    % Check if the maps L2[x,y]-->L2[x] and L2[x]-->L2[x,y] are defined by
    % full integral operators
    is_full_int(4) = false;
elseif ~all(all(isequal(cleanpoly(Pop.Ry2{2}-Pop.Ry2{3},Ctol),0))) ||...
        ~all(all(isequal(cleanpoly(Pop.R2y{2}-Pop.R2y{3},Ctol),0)))
    % Check if the maps L2[x,y]-->L2[y] and L2[y]-->L2[x,y] are defined by
    % full integral operators
    is_full_int(5);
end

% % Check that off-diagonal multipliers are zero.
is_no_mult = true(3,1);
if ~all(all(isequal(cleanpoly(Pop.Rx2{1},Ctol),0))) ||...
        ~all(all(isequal(cleanpoly(Pop.R2x{1},Ctol),0)))
    is_no_mult(1) = false;
elseif ~all(all(isequal(cleanpoly(Pop.Ry2{1},Ctol),0))) ||...
        ~all(all(isequal(cleanpoly(Pop.R2y{1},Ctol),0)))
    is_no_mult(2) = false;
elseif ~all(all(isequal(cleanpoly(Pop.R22{2,1},Ctol),0))) ||...
            ~all(all(isequal(cleanpoly(Pop.R22{3,1},Ctol),0))) ||...
                ~all(all(isequal(cleanpoly(Pop.R22{1,2},Ctol),0))) ||...
                    ~all(all(isequal(cleanpoly(Pop.R22{1,3},Ctol),0)))
    is_no_mult(3) = false;
end

% % Determine if the operator is separable
is_sep = all(is_full_int) && all(is_no_mult);

end
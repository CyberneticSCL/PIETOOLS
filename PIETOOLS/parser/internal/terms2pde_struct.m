function [PDE_out] = terms2pde_struct(terms_in)
% PDE_OUT = TERMS2PDE_STRUCT(TERMS_IN) converted an object of type 'terms'
% to an object of type 'pde_struct'. The 'terms' objects are built on the
% retired 'opvar' class, so the conversion is no longer supported, and the
% function throws an error.
%
% INPUT
% - terms_in:   'terms' class object;
%
% OUTPUT
% - PDE_out:    none, an error is thrown.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding DJ - 06/30/2024
% MMP, 11/19/2026: Release comment pass (comment-only).
% DJ, 11/25/2026: Throw an error, as 'terms' objects are built on the
%                   retired 'opvar' class; the conversion (DJ, 06/30/2024),
%                   and its subroutine is_LHS_term, are removed;

error("Objects of type 'terms' (a legacy input format) are no longer supported. Declare the PDE using 'pde_var' instead.")

end

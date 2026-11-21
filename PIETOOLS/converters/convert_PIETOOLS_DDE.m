%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_DDE.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DDF,PIE] = convert_PIETOOLS_DDE(DDE,out_type)
% This a setup routine for converting DDEs into DDFs and PIEs
% First converts DDE representation to minimal DDF representation: one
% delay channel r_i per delay tau(i), of the smallest dimension that
% reproduces the discrete and distributed delay terms of tau(i).
% Then converts DDF representation to PIE representation.
% If out_type=='ddf', only the DDF structure is computed and returned.
% If out_type=='pie', only PIE structure is returned (though DDF is also
% computed).
%
% A Partial Integral Equation is defined by 12 PI operators as
%
% Twop \dot{w}(t)+Tuop \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B1op w(t) + B2op u(t)
%                                           z(t)= C1op x(t) + D11op w(t) + D12op u(t)
%                                           y(t)= C2op x(t) + D21op w(t) + D22op u(t)
%
% The formulae were defined in (12) and (15)-(16) in the paper
% ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs''
% reference: https://arxiv.org/abs/1910.03881
%
% INPUT
%   DDE:      'struct', the DDE with K = numel(DDE.tau) delays
%     \dot x(t) = A0 x(t) + sum_i Ai{i} x(t-tau(i)) + sum_i int_{-tau(i)}^0 Adi{i}(s) x(t+s) ds
%               + B1 w(t) + (B1i, B1di likewise) + B2 u(t) + (B2i, B2di likewise)
%          z(t) = C1 x(t) + D11 w(t) + D12 u(t) + (C1i, C1di, D11i, D11di, D12i, D12di likewise)
%          y(t) = C2 x(t) + D21 w(t) + D22 u(t) + (C2i, C2di, D21i, D21di, D22i, D22di likewise)
%             fields tau (1 x K 'double', delays > 0), the 'double' matrices
%             A0, ..., D22, and 1 x K cells Ai, ..., D22i of 'double'
%             matrices and Adi, ..., D22di of 'double' or 'polynomial'
%             kernels; undefined fields are zero.
%   out_type: (optional) 'char', 'pie', 'ddf' or 'ddf_min' (the last two
%             equal); with out_type, at most one output.
%
% OUTPUT
%   DDF:  'struct', the minimal DDF (fields tau, A0, B1, B2, C1, C2, D11,
%         D12, D21, D22, Bv, D1v, D2v, and 1 x K cells Cr, Br1, Br2, Drv,
%         Cv, Cvd); the PIE instead when out_type = 'pie'.
%   PIE:  'pie_struct', the PIE of the DDF: 'copvar' operators on
%         R^nx x L2^nr[-1,0] in s1, nr the total dimension of the
%         delay channels.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding: DJ, 09/29/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.
%

% Check that desired conversion makes sense.
return_pie = false;
if nargin>=2
    if nargout>1
        error('At most one output is returned when an output system type is specified.')
    end
    if strcmpi(out_type,'pie')
        return_pie = true;
    elseif ~strcmpi(out_type,'ddf') && ~strcmpi(out_type,'ddf_min')
        error('Second argument must be one of ''pie'' or ''ddf''. PIETOOLS cannot convert DDEs to any other type.')
    end
end

% First convert to DDF.
DDF = minimize_PIETOOLS_DDE2DDF(DDE);
if nargout==1 && ~return_pie
    return
end

% Then convert DDF to PIE.
PIE = convert_PIETOOLS_DDF(DDF,'pie');
if return_pie
    % If the PIE is specifically requested, return the PIE as first
    % argument.
    DDF = PIE;
    return
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_NDS.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DDF,DDF_min,PIE] = convert_PIETOOLS_NDS(NDS,out_type)
% This a setup routine for converting NDSs into PIEs
% First converts NDS representation to DDF representation, with one
% delay channel r_i = [x; w; u; \dot x] per delay tau(i).
% Then minimizes the DDF representation: each channel is replaced by one
% of the smallest dimension with the same delayed terms.
% Then, converts DDF representation to PIE representation.
% If out_type=='ddf_max', the function will compute and return ONLY the
% non-minimized DDF representation.
% If out_type=='ddf' or out_type=='ddf_min', the function will also compute the
% minimized representation, and return ONLY this minimized representation.
% If out_type=='pie', the function will also compute the PIE representation,
% and return ONLY this PIE representation.
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
%   NDS:      'struct', the neutral delay system with K = numel(NDS.tau) delays
%     \dot x(t) = A0 x(t) + sum_i Ai{i} x(t-tau(i)) + sum_i int_{-tau(i)}^0 Adi{i}(s) x(t+s) ds
%               + sum_i Ei{i} \dot x(t-tau(i)) + sum_i int_{-tau(i)}^0 Edi{i}(s) \dot x(t+s) ds
%               + B1 w(t) + (B1i, B1di likewise) + B2 u(t) + (B2i, B2di likewise)
%          z(t) = C1 x(t) + D11 w(t) + D12 u(t) + (C1i, C1di, C1ei, C1dei, D11i, ..., D12di likewise)
%          y(t) = C2 x(t) + D21 w(t) + D22 u(t) + (C2i, C2di, C2ei, C2dei, D21i, ..., D22di likewise)
%             fields tau (1 x K 'double', delays > 0), the 'double' matrices
%             A0, ..., D22, and 1 x K cells of 'double' matrices (Ai, Ei,
%             C1ei, ...) and of 'double' or 'polynomial' kernels (Adi, Edi,
%             C1dei, ...); undefined fields are zero.
%   out_type: (optional) 'char', 'ddf_max', 'ddf', 'ddf_min' (equal to
%             'ddf') or 'pie'; with out_type, at most one output.
%
% OUTPUT
%   DDF:      'struct', the DDF of the NDS (fields tau, A0, B1, B2, C1, C2,
%             D11, D12, D21, D22, Bv, D1v, D2v, and 1 x K cells Cr, Br1,
%             Br2, Drv, Cv; no field Cvd); the minimal DDF or the PIE
%             instead when out_type is 'ddf', 'ddf_min' or 'pie'.
%   DDF_min:  'struct', the minimal DDF, the fields of DDF and the 1 x K
%             cell Cvd.
%   PIE:      'pie_struct', the PIE of the minimal DDF: 'copvar' operators
%             on R^nx x L2^nr[-1,0] in s1, nr the total dimension of the
%             delay channels.

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
% MMP, 11/19/2026: Release comment pass (comment-only).
%

% Check that desired conversion makes sense.
return_ddf = false;
return_pie = false;
if nargin>=2
    if nargout>1
        error('At most one output is returned when a system type for the output is specified.')
    end
    if strcmpi(out_type,'ddf') || strcmpi(out_type,'ddf_min')
        return_ddf = true;
    elseif strcmpi(out_type,'pie')
        return_pie = true;
    elseif ~strcmpi(out_type,'ddf_max')
        error('Second argument must be one of ''pie'',''ddf'', or ''ddf_max''. PIETOOLS cannot convert NDSs to any other type.')
    end
end

% First convert to DDF.
DDF = convert_PIETOOLS_NDS2DDF(NDS);
if nargout==1 && ~return_ddf && ~return_pie
    return
end

% Then minimize the DDF.
DDF_min = minimize_PIETOOLS_DDF(DDF);
if nargout==2
    return
elseif return_ddf
    % If a minimal DDF representation is explicitly requested, return this
    % as first and only argument.
    DDF = DDF_min;
    return
end

% Then convert DDF to PIE.
PIE = convert_PIETOOLS_DDF(DDF_min,'pie');
if return_pie
    % If a PIE representation is explicitly requested, return this as first
    % and only argument.
    DDF = PIE;
    return
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_NDS.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DDF,DDF_min,PIE] = convert_PIETOOLS_NDS(NDS,out_type)
% This a setup routine for converting NDSs into PIEs
% First converts NDS representation to DDF representation using
% "convert_PIETOOLS_NDS2DDF".
% Then minimizes the DDF representation using
% "minimize_PIETOOLS_DDF".
% Then, converts DDF representation to PIE representation using
% "convert_PIETOOLS_DDF".
% If type=='ddf_max', the function will compute and return ONLY the 
% non-minimized DDF representation.
% If type=='ddf' or type=='ddf_min', the function will also compute the
% minimized representation, and return ONLY this minimized representation.
% If type=='pie', the function will also compute the PIE representation,
% and return ONLY this PIE representation.
% 
% A Partial Integral Equation is defined by 12 PI operators as
%
% Twop \dot{w}(t)+Tuop \dot{u}(t)+Top \dot{x}(t)= Aop x(t) + B1op u(t)+ + B2op w(t)
%                                           z(t)= C1op x(t) + D11op u(t)+ + D12op w(t)
%                                           y(t)= C2op x(t) + D21op u(t)+ + D22op w(t)
%
% The formulae were defined in (12) and (15)-(16) in the paper
% ''Representation of Networks and Systems with Delay: DDEs, DDFs, ODE-PDEs and PIEs''
% reference: https://arxiv.org/abs/1910.03881

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
% Initial coding: DJ, 09/29/2022
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
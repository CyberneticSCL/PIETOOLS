%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_DDE.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DDF,PIE] = convert_PIETOOLS_DDE(DDE,out_type)
% This a setup routine for converting DDEs into DDFs and PIEs
% First converts DDE representation to minimal DDF representation using
% "minimize_PIETOOLS_DDE2DDF".
% Then converts DDF representation to PIE representation using
% "convert_PIETOOLS_DDF".
% If type=='ddf', only the DDF structure is computed and returned.
% If type=='pie', only PIE structure is returned (though DDF is also
% computed).
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
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding: DJ, 09/29/2022
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
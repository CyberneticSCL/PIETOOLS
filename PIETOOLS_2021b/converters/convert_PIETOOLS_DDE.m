%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_DDE.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PIE,DDF] = convert_PIETOOLS_DDE(DDE,out_type)
% This a setup routine for converting DDEs into DDFs and PIEs
% First converts DDE representation to minimal DDF representation using
% "minimize_PIETOOLS_DDE2DDF".
% If type=='ddf', the DDF structure will be returned.
% If type=='pie', converts DDF representation to PIE representation using
% "convert_PIETOOLS_DDF".
% The function returns the PIE representation, and the intermediate
% DDF representation.
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
return_ddf = false;
if nargin>=2 && (strcmpi(out_type,'ddf') || strcmpi(out_type,'ddf_min'))
    return_ddf = true;
    if nargout>1
        error('At most 1 output is supported for DDE to DDF conversion.')
    end
elseif nargin>=2 && ~strcmpi(out_type,'pie')
    error('Second argument must be one of ''pie'' or ''ddf''. PIETOOLS cannot convert DDEs to any other type.')
end

% First convert to DDF.
DDF = minimize_PIETOOLS_DDE2DDF(DDE);
if return_ddf
    PIE = DDF;
    return
end

% Then convert DDF to PIE.
PIE = convert_PIETOOLS_DDF(DDF);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_NDS.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PIE,DDF_min,DDF] = convert_PIETOOLS_NDS(NDS,out_type)
% This a setup routine for converting NDSs into PIEs
% First converts NDS representation to DDF representation using
% "convert_PIETOOLS_NDS2DDF".
% If type=='ddf_max', the function will return the DDF representation as
% first and only argument.
% If type=='ddf' or type=='PIE', minimizes the DDF representation using
% "minimize_PIETOOLS_DDF".
% If type=='ddf', the function will return the minimized DDF representation
% as first argument, and the original DDF representation as second
% argument.
% If type=='pie', converts DDF representation to PIE representation using
% "convert_PIETOOLS_DDF".
% The function returns the PIE representation, the minimized DDF
% representation, and the non-minimized DDF representation.
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
return_ddf_max = false;
if nargin>=2 && strcmpi(out_type,'ddf_max')
    return_ddf_max = true;
    if nargout>2
        error('At most 1 output2 is supported for NDS to non-minimal DDF conversion.')
    end
elseif nargin>=2 && (strcmpi(out_type,'ddf') || strcmpi(out_type,'ddf_min'))
    return_ddf = true;
    if nargout>2
        error('At most 2 outputs are supported for NDS to DDF conversion.')
    end
elseif nargin>=2 && ~strcmpi(out_type,'pie')
    error('Second argument must be one of ''pie'',''ddf'', or ''ddf_max''. PIETOOLS cannot convert NDs to any other type.')
end

% First convert to DDF.
DDF = convert_PIETOOLS_NDS2DDF(NDS);
if return_ddf_max
    PIE = DDF;
    return
end

% Then minimize the DDF.
DDF_min = minimize_PIETOOLS_DDF(DDF);
if return_ddf
    PIE = DDF_min;
    DDF_min = DDF;
    return
end

% Then convert DDF to PIE.
PIE = convert_PIETOOLS_DDF(DDF_min);

end
function [prog, varargout] = lpisolve(PIE,settings,lpi,fn_hnd)
if nargin<3
    error('Insufficient number of arguments. Correct syntax: lpisolve(PIE,settings,lpi,function_handle)');
end
if isempty(settings)
    settings = lpisettings('light');
end
if ~ismember(lpi,{'stability','stability-dual','l2-gain','l2-gain-dual','hinf-observer','hinf-controller','custom'})
    error("Unknown lpi type. Use 'custom' option and pass the function handle for the lpi to be solved");
end
if strcmp(lpi,'custom') && nargin~=4
    error("lpi type was set to 'custom', however, no function handle was passed. Correct syntax: lpisolve(PIE,settings,lpi,function_handle)");
end



switch lpi
    case 'stability'
        [prog, P] = PIETOOLS_stability(PIE,settings);
        varargout{1} = P;
    case 'stability-dual'
        [prog, P] = PIETOOLS_stability_dual(PIE,settings);
        varargout{1} = P;
    case 'l2-gain'
        [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);
        varargout{1} = P; varargout{2} = gamma;
    case 'l2-gain-dual'
        [prog, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,settings);
        varargout{1} = P; varargout{2} = gamma;
    case 'hinf-observer'
        [prog, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);
        varargout{1} = P; varargout{2} = L; varargout{3} = gamma; varargout{4} = P; varargout{5} = Z;
    case 'hinf-controller'
        [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);
        varargout{1} = P; varargout{2} = K; varargout{3} = gamma; varargout{4} = P; varargout{5} = Z;
    case 'custom'
        [prog, vout] = fn_hnd(PIE,settings);
        varargout = cell(1,length(vout));
        for i=1:length(vout)
            varargout{i} = vout{i};
        end
end
end
function [prog, varargout] = lpisolve(PIE,settings,lpi)
if nargin<3
    error('Insufficient number of arguments. Correct syntax: lpisolve(PIE,settings,lpi)');
end
if isempty(settings)
    settings = lpisettings('light');
elseif ~isa(settings,'string')&&~isa(settings,'struct')
    error("Settings must either be a string value or a settings structure similar to the output of lpisettings() function.")
end
if ~ismember(lpi,{'stability','stability-dual','l2-gain','l2-gain-dual','hinf-observer','hinf-controller','custom'})
    error("Unknown lpi type. Use 'custom' option and pass the function handle for the lpi to be solved");
elseif ~isa(lpi,'function_handle')
    error("lpi must be a string value or a function handle.");
end

if isa(settings,'string')
    settings = lpisettings(settings);
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
    otherwise
        [prog, vout] = lpi(PIE,settings);
        varargout = cell(1,length(vout));
        for i=1:length(vout)
            varargout{i} = vout{i};
        end
end
end
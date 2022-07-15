function sttgns = lpisettings(type,dim,derivative_strictness,simplify)
arguments
    type {mustBeMember(type,{'light','heavy','veryheavy','stripped','extreme','custom'})}='light';
    dim {mustBeInteger,mustBeInRange(dim,1,2)} = 1;
    derivative_strictness {mustBeNonnegative}= 0;
    simplify {mustBeMember(simplify,{'','psimplify'})}= '';
end

if dim==1
    switch type
        case 'light'
            sttgns = settings_PIETOOLS_light;
        case 'custom'
            sttgns = settings_PIETOOLS_custom;
        case 'extreme'
            sttgns = settings_PIETOOLS_extreme;            
        case 'heavy'
            sttgns = settings_PIETOOLS_heavy;
        case 'stripped'
            sttgns = settings_PIETOOLS_stripped;
        case 'veryheavy'
            sttgns = settings_PIETOOLS_veryheavy;
        otherwise
            warning('Unknown settings parameter requested. Defaulting to light settings');
    end
else
    switch type
    end
end

if strcmp(simplify,'psimplify')
    sttgns.sos_opts.simplify = 1; % Use psimplify
else
    sttgns.sos_opts.simplify = 0; % Dont Use psimplify
end

sttgns.eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
sttgns.eppos2 = 1*1e-6;   % Positivity of Lyapunov Function with respect to spatially distributed states

if derivative_strictness
    sttgns.epneg = derivative_strictness;    % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
end
end
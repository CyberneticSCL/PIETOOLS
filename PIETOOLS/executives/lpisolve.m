function [prog, varargout] = lpisolve(PIE,settings,lpi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a wrapper function for various executives in PIETOOLS. 
% Input -
% PIE: pie_struct class object
% settings: lpisettings() object or string specifying settings type 'light','heavy','veryheavy','stripped','extreme','custom'
% lpi: lpi executive function handle or string specifying lpi type, 'stability', 'stability-dual','l2gain','l2gain-dual', 'hinf-observer','hinf-controller'
% Output -
% prog: solved sosprogram() structure
% varargout: output variables from different lpi types
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if nargin<3
    error('Insufficient number of arguments. Correct syntax: lpisolve(PIE,settings,lpi)');
end
if isempty(settings)
    settings = lpisettings('light');
elseif ~isa(settings,'string')&&~isa(settings,'struct')
    error("Settings must either be a string value or a settings structure similar to the output of lpisettings() function.")
end
if (isa(lpi,'string')||isa(lpi,'char'))&&ismember(lpi,{'stability','stability-dual','l2gain','l2gain-dual','hinf-observer','hinf-controller','custom'})
    % do nothing
elseif ~isa(lpi,'function_handle')
    error("lpi must be a string value or a function handle.");
else
    error("Unknown lpi type");
end

if isa(settings,'string')
    settings = lpisettings(settings);
end

if isa(PIE,'sys')
    PIE = PIE.params;
end

switch lpi
    case 'stability'
        [prog, P] = PIETOOLS_stability(PIE,settings);
        varargout{1} = P;
    case 'stability-dual'
        [prog, P] = PIETOOLS_stability_dual(PIE,settings);
        varargout{1} = P;
    case 'l2gain'
        [prog, P, gamma] = PIETOOLS_Hinf_gain(PIE,settings);
        varargout{1} = P; varargout{2} = gamma;
    case 'l2gain-dual'
        [prog, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,settings);
        varargout{1} = P; varargout{2} = gamma;
    case 'hinf-observer'
        [prog, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,settings);
        varargout{1} = L; varargout{2} = gamma; varargout{3} = P; varargout{4} = Z; 
    case 'hinf-controller'
        [prog, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,settings);
        varargout{1} = K; varargout{2} = gamma; varargout{3} = P; varargout{4} = Z; 
    otherwise
        [prog, vout] = lpi(PIE,settings);
        varargout = cell(1,length(vout));
        for i=1:length(vout)
            varargout{i} = vout{i};
        end
end
end
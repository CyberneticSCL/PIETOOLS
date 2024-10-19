function [prog_sol, varargout] = lpisolve(lpi,opts,PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a wrapper function for various executives in PIETOOLS. 
% INPUT
% - lpi:        Either an LPI program structure (based on SOSTOOLS
%               structure) specifying an LPI, or a 'char' or 'struct'
%               specifying one of the executives to call, namely one of:
%               'stability', 'stability-dual','l2gain','l2gain-dual',
%                   'hinf-observer','hinf-controller';
% - opts:       Either options to be passed to SOSTOOLS for solving the
%               optimization program specified by lpi, or if 'lpi'
%               corresponds to one of the predefined executives, a struct
%               specifying settings to use in declaring and solving that
%               LPI. In the latter case, can also be a 'char' or 'struct'
%               specifying desired settings to use, one of:
%               'light','heavy','veryheavy','stripped','extreme','custom';
% - PIE:        Optional argument in case a pre-defined LPI is to be solved
%               specifying a desired PIE for which to solve that LPI, 
%               specified as pie_struct class object.
%
% OUTPUT
% - prog:       solved sosprogram() structure for desired LPI
% - varargout:  Additional outputs returned from pre-defined LPI
%                   executives.
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
%
% DJ, 10/18/2024 - Initial coding


% % % There are two ways to use this function:
if isa(lpi,'struct') && isfield(lpi,'vartable')
    %% OPTION 1: lpi is an LPI program structure
    
    % % In this case, 'settings' should correspond to SOSTOOLS settings for
    % % 'sos_solve', and 'PIE' should be empty
    if nargin>=3
        error('Too many input arguments.')
    end
    % % We output only the solved optimization program structure.
    if nargout>1
        error('Too many output arguments.')
    end
    % % Solve using 'sossolve' with desired settings, if specified.
    if nargin>=2
        prog_sol = sossolve(lpi,opts);
    else
        prog_sol = sossolve(lpi);
    end
else
    %% OPTION 2: lpi is one of the pre-built LPI executives.
    
    if nargin<3
        error('Insufficient number of arguments. Correct syntax: lpisolve(PIE,settings,lpi)');
    end
    
    % Check which LPI to solve
    if isa(lpi,'function_handle')
        % Assume the input 'lpi' refers to a function that can be called as
        %   [prog, vout] = lpi(PIE,settings);
        % Do nothing here...
    elseif ~(isa(lpi,'string') || isa(lpi,'char'))
        error("The input 'lpi' must be a string value or a function handle.");
    elseif ~ismember(lpi,{'stability','stability-dual','l2gain','l2gain-dual','hinf-observer','hinf-controller','custom'})
        error("Specified LPI type must be one of 'stability', 'stability-dual', 'l2gain', 'l2gain-dual', 'hinf-observer', or 'hinf-controller'.")
    end
    
    % Check what settings to use
    if isempty(opts)
        opts = lpisettings('light');
    elseif ~(isa(opts,'string')||isa(opts,'char'))&&~isa(opts,'struct')
        error("Settings must either be a string value or a settings structure similar to the output of lpisettings() function.")
    end
    if isa(opts,'string')||isa(opts,'char')
        opts = lpisettings(opts);
    end
    
    % Check for what system to solve
    if isa(PIE,'sys') && strcmpi(PIE.type,'pie')
        PIE = PIE.params;
    elseif (isa(PIE,'sys') && strcmpi(PIE.type,'pde')) || isa(PIE,'pde_struct')
        % Convert PDE to PIE
        PIE = convert(PIE,'pie');
    end
    
    % Run the desired executive for the specified system.
    switch lpi
        case 'stability'
            [prog_sol, P] = PIETOOLS_stability(PIE,opts);
            varargout{1} = P;
        case 'stability-dual'
            [prog_sol, P] = PIETOOLS_stability_dual(PIE,opts);
            varargout{1} = P;
        case 'l2gain'
            [prog_sol, P, gamma] = PIETOOLS_Hinf_gain(PIE,opts);
            varargout{1} = P; varargout{2} = gamma;
        case 'l2gain-dual'
            [prog_sol, P, gamma] = PIETOOLS_Hinf_gain_dual(PIE,opts);
            varargout{1} = P; varargout{2} = gamma;
        case 'hinf-observer'
            [prog_sol, L, gamma, P, Z] = PIETOOLS_Hinf_estimator(PIE,opts);
            varargout{1} = L; varargout{2} = gamma; varargout{3} = P; varargout{4} = Z; 
        case 'hinf-controller'
            [prog_sol, K, gamma, P, Z] = PIETOOLS_Hinf_control(PIE,opts);
            varargout{1} = K; varargout{2} = gamma; varargout{3} = P; varargout{4} = Z; 
        otherwise
            [prog_sol, vout] = lpi(PIE,opts);
            varargout = cell(1,length(vout));
            for i=1:length(vout)
                varargout{i} = vout{i};
            end
    end
end


end
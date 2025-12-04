function [prog_sol, varargout] = lpiscript(PIE,lpi,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a wrapper function for various executives in PIETOOLS. 
% INPUT
% - PIE:        A 'pie_struct' object specifying the PIE system for which
%               to solve the desired LPI;
% - lpi:        A 'char' or 'string' specifying one of the executives to 
%               call, namely one of:
%               'stability', 'stability-dual','l2gain','l2gain-dual',
%                 'h2norm','h2norm-dual','hinf-observer','hinf-controller',
%                   'h2-observer','h2-controller';
% - opts:       A 'struct' specifying settings to use in declaring and 
%               solving the desired pre-defined LPI. Alternatively,
%               this field can also be a 'char' or 'string' specifying
%               desired settings to use, namely one of:
%               'extreme','stripped','light','heavy','veryheavy','custom';
%
% OUTPUT
% - prog:       solved LPI program structure for desired LPI;
% - varargout:  Additional outputs returned from pre-defined LPI
%                   executives;
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024 PIETOOLS Team
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
% DJ, 11/01/2024: Initial coding;
% DJ, 01/06/2025: Add support for call to H2 executives; 
% DB 01/07/2025: Fix output of H2 executives

% % % There are two ways to use this function: 
if nargin<3
    error('Insufficient number of arguments. Correct syntax: lpisolve(PIE,lpi,settings)');
end

% Check which LPI to solve
if isa(lpi,'function_handle')
    % Assume the input 'lpi' refers to a function that can be called as
    %   [prog, vout] = lpi(PIE,settings);
    % Do nothing here...
elseif ~(isa(lpi,'string') || isa(lpi,'char'))
    error("The input 'lpi' must be a string value or a function handle.");
elseif ~ismember(lpi,{'stability','stability-dual',...
                        'l2gain','l2gain-dual','h2norm','h2norm-dual',...
                        'hinf-observer','hinf-controller','h2-observer','h2-controller','custom'})
    error("Specified LPI type must be one of 'stability', 'stability-dual'," + ...
            " 'l2gain', 'l2gain-dual', 'h2norm', 'h2norm-dual',"+...
             " 'hinf-observer', 'hinf-controller', 'h2-observer', or 'h2-controller'.")
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
elseif ~isa(PIE,'pie_struct')
    error('PIE system for which to run the executive should be specified as object of type ''pie_struct''.')
end

% Run the desired executive for the specified system.
switch lpi
    case 'stability'
        [prog_sol, P] = PIETOOLS_stability_PIE2PDE(PIE,opts);
        varargout{1} = P;
    case 'stability-dual'
        [prog_sol, P] = PIETOOLS_stability_dual(PIE,opts);
        varargout{1} = P;
    case 'l2gain'
        [prog_sol, P, gam] = PIETOOLS_Hinf_gain(PIE,opts);
        varargout{1} = P; varargout{2} = gam;
    case 'l2gain-dual'
        [prog_sol, P, gam] = PIETOOLS_Hinf_gain_dual(PIE,opts);
        varargout{1} = P; varargout{2} = gam;
    case 'h2norm'
        [prog_sol, W, gam, R, Q]= PIETOOLS_H2_norm_o(PIE,opts);
        varargout{1} = W; varargout{2} = gam; varargout{3} = R; varargout{4} = Q;
    case 'h2norm-dual'
        [prog_sol, W, gam, R, Q] = PIETOOLS_H2_norm_c(PIE,opts);
        varargout{1} = W; varargout{2} = gam; varargout{3} = R; varargout{4} = Q;
    case 'hinf-observer'
        [prog_sol, L, gam, P, Z] = PIETOOLS_Hinf_estimator(PIE,opts);
        varargout{1} = L; varargout{2} = gam; varargout{3} = P; varargout{4} = Z; 
    case 'hinf-controller'
        [prog_sol, K, gam, P, Z] = PIETOOLS_Hinf_control(PIE,opts);
        varargout{1} = K; varargout{2} = gam; varargout{3} = P; varargout{4} = Z; 
    case 'h2-observer'
        [prog_sol, L, gam, P, Z,W] = PIETOOLS_H2_estimator(PIE,opts);
        varargout{1} = L; varargout{2} = gam; varargout{3} = P; varargout{4} = Z;varargout{5} = W;
    case 'h2-controller'
        [prog_sol, K, gam, P, Z,W] = PIETOOLS_H2_control(PIE,opts);
        varargout{1} = K;  varargout{2} = gam; varargout{3} = P; varargout{4} = Z;varargout{5} = W;
    otherwise
        [prog_sol, vout] = lpi(PIE,opts);
        varargout = cell(1,length(vout));
        for i=1:length(vout)
            varargout{i} = vout{i};
        end
end


end
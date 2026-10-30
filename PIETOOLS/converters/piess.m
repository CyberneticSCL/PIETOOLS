function PIE = piess(Top,Aop,Bop,Cop,Dop)
% PIE = PIESS(TOP,AOP,BOP,COP,DOP) takes a set of PI operators, and
% constructs the associated PIE representation. Declaring Top through Dop
% as 'copvar' objects, the output PIE is a 'pie_struct' object
% representing a "5-PIE",
%
%   d/dt Top*v(t) = Aop*v(t) + Bop*w(t)
%            z(t) = Cop*v(t) + Dop*w(t)
%
% where v(t) is the fundamental state, w(t) the exogenous input
% (disturbance), and z(t) the regulated output, all possibly vector-valued
% and involving coupled finite- and infinite-dimensional components,
% determined by the dimensions of the specified operators.
% Alternatively, if cells of operators are specified as
%
%       Top = {Tvop,Twop,Tuop},     Bop = {Bwop,Buop};
%       Cop = {Czop;Cuop},          Dop = {Dzwop,Dzuop; Dywop,Dyup};
% 
% the output PIE is a 'pie_struct' object representing the "12-PIE"
% 
%   d/dt (Top*v(t)+Twop*w(t)+Tuop*u(t)) = Aop*v(t) +Bwop*w(t) +Buop*u(t)
%                                  z(t) = Czop*v(t)+Dzwop*w(t)+Dzuop*u(t)
%                                  y(t) = Cyop*v(t)+Dywop*w(t)+Dyuop*u(t)
% 
% where now u(t) is the controlled input, and y(t) the observed output,
% with dimensions and domain again determined by the input operators.
%
% INPUT
% - Top:    nv x nv 'copvar' object representing the
%           PI operator defining the left-hand side of a 5-PIE.
%           Alternatively, a 1x3 cell {Tvop,Twop,Tuop} with all elements 
%           'copvar' objects of dimensions 
%           {nv x nv, nv x nw, nv x nu}, representing the left-hand
%           side of a 12-PI. The operators Twop and Tuop can also be set
%           equal to [], in which case they will default to 0 operators
%           of suitable dimensions;
% - Aop:    nv x nv 'copvar' object representing the state
%           operator in both 5-PIE and 12-PIE dynamics.
% - Bop:    nv x nw 'copvar' object representing the
%           input PI operator defining the contribution of disturbances w 
%           to the PDE dynamics.
%           Alternatively, a 1x2 cell {Bwop,Buop} with both elements 
%           'copvar' objects of dimensions 
%           {nv x nw, nv x nu}, representing the contribution of 
%           disturbances and controlled inputs to the PDE dynamics;
%           Defaults to [], returning 0 operator(s) of suitable dimensions;
% - Cop:    nz x nv 'copvar' object representing the
%           output PI operator defining the contribution of the state v 
%           to the regulated output z.
%           Alternatively, a 2x1 cell {Czop;Cyop} with both elements 
%           'copvar' objects of dimensions 
%           {nz x nv, ny x nv}, representing the contribution of 
%           the state to the regulated and controlled output.
%           Defaults to [], returning 0 operator(s) of suitable dimensions;
% - Dop:    nz x nw 'copvar' object representing the
%           feedthrough PI operator defining the effect of the 
%           disturbance w to the regulated output z.
%           Alternatively, a 2x2 cell {Dzwop, Dzuop; Dywop, Dyuop} with 
%           all elements 'copvar' objects of dimensions 
%           {nz x nw, nz x nu; ny x nw, ny x nu}, representing the
%           feedthrough from each type of input to each type of output.
%           Defaults to [], returning 0 operator(s) of suitable dimensions;
%
% OUTPUT
% - PIE:    'pie_struct' object representing the (linear) PIE defined by
%           the specified operators. The specified operators will be stored
%           in the PIE structure under the appropriate fields, so that e.g.
%           PIE.T = Top, PIE.Tw = Twop, PIE.A = Aop, etc.. The PIE is
%           initialized, so that the operators are 'copvar' objects.
%
% NOTES
% - The operators are 'copvar' objects, the type of the operators of a
% 'pie_struct' (e.g. as returned by 'convert'). If instead at least one is
% an 'opvar' or 'opvar2d' object, and none is a 'copvar', the PIE is
% constructed by 'piess_opvar'. Otherwise, 'opvar' and 'opvar2d' objects
% are converted to 'copvar' objects by INITIALIZE.
% - Matrices are interpreted as constant multipliers between the spaces of
% the 'copvar' operators acting on the same components, or between R^n if
% there are none. Top and Aop must be specified, unless the operators are
% handled by 'piess_opvar', which sets any unspecified operator to 0.
% - Only linear PIEs, of type 'pie_struct', are supported.
% - Passing a single argument "PIE = piess(Top)", the input Top will be
% interpreted as operator Dzwop from disturbance to regulated output, and a
% PIE will be constructed as just z(t) = Dzwop*w(t). However, PIETOOLS
% currently offers no functionality for actually working with such systems.
% - PIEs with algebraic constraints (K, Kw, Ku) are currently not
% supported.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piess
%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 12/20/2024: Initial coding;
% DJ, 12/30/2024: Suppress feedthrough warning if Dop is empty;
% DJ, 10/29/2026: Update to new pie_struct structure;
% DJ (CC), 10/30/2026: Build the PIE from 'copvar' operators, the type of
%                   the operators of a pie_struct. The construction from
%                   'opvar' and 'opvar2d' operators moved to piess_opvar;


% % % Process the inputs
% % Check how may input arguments are passed.
if nargin==0
    % Return empty pie_struct structure.
    PIE = pie_struct();
    return
elseif nargin==1
    % Just as Matlab's `ss`, interpret single argument as direct map from
    % input to output.
    Dop = Top;
    Top = [];   Aop = [];   Bop = [];   Cop = [];
elseif nargin==2
    % Interpret two arguments as autonomous system d/dt(Top*v) = Aop*v.
    Bop = [];   Cop = [];   Dop = [];
elseif nargin==3
    % Interpret three arguments as system with disturbance,
    %   d/dt(Top*v) = Aop*v + Bop*w
    Cop = [];   Dop = [];
elseif nargin==4
    % Interpret three arguments as system with disturbance and output,
    %   d/dt(Top*v) = Aop*v + Bop*w
    %             z = Cop*v
    Dop = [];
elseif nargin>5
    error('Too many input arguments.')
end

% % Extract operators defining the PIE.
% Types of object accepted as operator ('copvar' is the default).           % DJ (CC), 10/30/2026
is_op = @(P) isa(P,'copvar') || isa(P,'opvar') || isa(P,'opvar2d') || isa(P,'polynomial') || isa(P,'double');   % DJ (CC), 10/30/2026
% Operators Top s.t u = Top*v + Twop*w + Tuop*u
Twop = [];      Tuop = [];
if isa(Top,'cell')
    if isscalar(Top)
        Top = Top{1};
    elseif numel(Top)==2
        Twop = Top{2};
        Top = Top{1};
    elseif numel(Top)==3
        Tuop = Top{3};
        Twop = Top{2};
        Top = Top{1};
    else
        error("Operators Top defining the fundamental map should be defined as a single 'opvar' object or 1x3 cell {Top,Twop,Tuop} of 'opvar' objects.")
    end
end
if isequal(Top,0) || isempty(Top)
    Top = [];
%   elseif ~isa(Top,'opvar') && ~isa(Top,'opvar2d') && ~isa(Top,'polynomial') && ~isa(Top,'double')  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Top)                                                          % DJ (CC), 10/30/2026
    error("Operators Top defining the fundamental map should be defined as a single 'opvar' object or 1x3 cell {Top,Twop,Tuop} of 'opvar' objects.")
end
if (isa(Twop,'double') && isequal(Twop,0)) || isempty(Twop)
    Twop = [];
%   elseif ~isa(Twop,'opvar') && ~isa(Twop,'opvar2d') && ~isa(Twop,'polynomial') && ~isa(Twop,'double')  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Twop)                                                         % DJ (CC), 10/30/2026
    error("Operators Top defining the fundamental map should be defined as a single 'opvar' object or 1x3 cell {Top,Twop,Tuop} of 'opvar' objects.")
end
if (isa(Tuop,'double') && isequal(Tuop,0)) || isempty(Tuop)
    Tuop = [];
%   elseif ~isa(Tuop,'opvar') && ~isa(Tuop,'opvar2d') && ~isa(Tuop,'polynomial') && ~isa(Tuop,'double')  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Tuop)                                                         % DJ (CC), 10/30/2026
    error("Operators Top defining the fundamental map should be defined as a single 'opvar' object or 1x3 cell {Top,Twop,Tuop} of 'opvar' objects.")
end
% Operator Aop
if isa(Aop,'cell')
    if ~isscalar(Aop)
        error("Operator Aop defining the PIE dynamics should be defined as a single 'opvar' object.")
    else
        Aop = Aop{1};
    end
end
if (isa(Aop,'double') && isequal(Aop,0)) || isempty(Aop)
    Aop = [];
%   elseif ~isa(Aop,'opvar') && ~isa(Aop,'opvar2d') && ~isa(Aop,'polynomial') && ~isa(Aop,'double')  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Aop)                                                          % DJ (CC), 10/30/2026
    error("Operator Aop defining the PIE dynamics should be defined as a single 'opvar' object.")
end
% Operator Bop
Bwop = Bop;     Buop = [];      % assume no controlled input by default
if isa(Bop,'cell')
    if isscalar(Bop)
        Bwop = Bop{1};
    elseif numel(Bop)==2
        Bwop = Bop{1};
        Buop = Bop{2};
    else
        error("Input operators Bop should be defined as a single 'opvar' object Bwop or 1x2 cell {Bwop,Buop} of 'opvar' objects.")
    end
end
if (isa(Bwop,'double') && isequal(Bwop,0)) || isempty(Bwop)
    Bwop = [];
%   elseif ~isa(Bwop,'opvar') && ~isa(Bwop,'opvar2d') && ~isa(Bwop,'polynomial') && ~isa(Bwop,'double')  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Bwop)                                                         % DJ (CC), 10/30/2026
    error("Input operators Bop should be defined as a single 'opvar' object Bwop or 1x2 cell {Bwop,Buop} of 'opvar' objects.")
end
if (isa(Buop,'double') && isequal(Buop,0)) || isempty(Buop)
    Buop = [];
%   elseif ~isa(Buop,'opvar') && ~isa(Buop,'opvar2d') && ~isa(Buop,'polynomial') && ~isa(Buop,'double')  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Buop)                                                         % DJ (CC), 10/30/2026
    error("Input operators Bop should be defined as a single 'opvar' object Bwop or 1x2 cell {Bwop,Buop} of 'opvar' objects.")
end
% Operator Cop
Czop = Cop;     Cyop = [];      % assume no controlled input by default
if isa(Cop,'cell')
    if isscalar(Cop)
        Czop = Cop{1};
    elseif numel(Cop)==2
        Czop = Cop{1};
        Cyop = Cop{2};
    else
        error("Output operators Cop should be defined as a single 'opvar' object Czop or 2x1 cell {Czop;Cyop} of 'opvar' objects.")
    end
end
if (isa(Czop,'double') && isequal(Czop,0)) || isempty(Czop)
    Czop = [];
%   elseif ~isa(Czop,'opvar') && ~isa(Czop,'opvar2d') && ~isa(Czop,'polynomial') && ~isa(Czop,'double')  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Czop)                                                         % DJ (CC), 10/30/2026
    error("Output operators Cop should be defined as a single 'opvar' object Czop or 2x1 cell {Czop;Cyop} of 'opvar' objects.")
end
if (isa(Cyop,'double') && isequal(Cyop,0)) || isempty(Cyop)
    Cyop = [];
%   elseif ~isa(Cyop,'opvar') && ~isa(Cyop,'opvar2d') && ~isa(Cyop,'polynomial') && ~isa(Cyop,'double')  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Cyop)                                                         % DJ (CC), 10/30/2026
    error("Output operators Cop should be defined as a single 'opvar' object Czop or 2x1 cell {Czop;Cyop} of 'opvar' objects.")
end
% Operator Dop
Dzwop = Dop;    Dzuop = [];
Dywop = [];     Dyuop = [];
if isa(Dop,'cell')
    if isscalar(Dop)
        Dzwop = Dop{1};
    elseif all(size(Dop)==[1,2])
        Dzwop = Dop{1};     Dzuop = Dop{2};
    elseif all(size(Dop)==[2,1])
        Dzwop = Dop{1};     Dywop = Dop{2};
    elseif all(size(Dop)==[2,2])
        Dzwop = Dop{1,1};   Dzuop = Dop{1,2};
        Dywop = Dop{2,1};   Dyuop = Dop{2,2};
    else
        error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
    end
%   elseif isa(Dop,'opvar') || isa(Dop,'opvar2d') || isa(Dop,'double') || isa(Dop,'polynomial')  % DJ (CC), 10/30/2026 (was)
elseif is_op(Dop)                                                           % DJ (CC), 10/30/2026
    % If a single operator is specified, check if only a single input or
    % output operator is specified as well.
    if isempty(Bwop) && ~isempty(Buop)
        if isempty(Czop) && ~isempty(Cyop)
            Dyuop = Dop;    Dzwop = [];
        else
            if ~isempty(Cyop) && ~isempty(Dop)                              % DJ, 12/30/2024
                warning("Assuming feedthrough operator maps to regulated output y.")
            end
            Dzuop = Dop;    Dzwop = [];
        end
    else 
        if ~isempty(Buop) && ~isempty(Dop)                                  % DJ, 12/30/2024
            warning("Assuming feedthrough operator maps disturbance w.")
        end
        if isempty(Czop) && ~isempty(Cyop)
            Dywop = Dop;    Dzwop = [];
        elseif ~isempty(Cyop) && ~isempty(Dop)                              % DJ, 12/30/2024
            warning("Assuming feedthrough operator maps to regulated output y.")
        end
    end
else
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end
if isa(Dzwop,'double') && (isempty(Dzwop) || isequal(Dzwop,0))
    Dzwop = [];
%   elseif ~(isa(Dzwop,'opvar') || isa(Dzwop,'opvar2d') || isa(Dzwop,'double') || isa(Dzwop,'polynomial'))  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Dzwop)                                                        % DJ (CC), 10/30/2026
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end
if isa(Dzuop,'double') && (isempty(Dzuop) || isequal(Dzuop,0))
    Dzuop = [];
%   elseif ~(isa(Dzuop,'opvar') || isa(Dzuop,'opvar2d') || isa(Dzuop,'double') || isa(Dzuop,'polynomial'))  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Dzuop)                                                        % DJ (CC), 10/30/2026
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end
if isa(Dywop,'double') && (isempty(Dywop) || isequal(Dywop,0))
    Dywop = [];
%   elseif ~(isa(Dywop,'opvar') || isa(Dywop,'opvar2d') || isa(Dywop,'double') || isa(Dywop,'polynomial'))  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Dywop)                                                        % DJ (CC), 10/30/2026
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end
if isa(Dyuop,'double') && (isempty(Dyuop) ||  isequal(Dyuop,0))
    Dyuop = [];
%   elseif ~(isa(Dyuop,'opvar') || isa(Dyuop,'opvar2d') || isa(Dyuop,'double') || isa(Dyuop,'polynomial'))  % DJ (CC), 10/30/2026 (was)
elseif ~is_op(Dyuop)                                                        % DJ (CC), 10/30/2026
    error("Feedthrough operators Dop should be declared as a single 'opvar' object Dzwop or 2x2 cell {Dzwop,Dzuop;Dywop,Dyuop} of 'opvar' objects.")
end

% % % opvar operators, and no copvar: piess_opvar.                          % DJ (CC), 10/30/2026
Pops = {Top,Twop,Tuop; Aop,Bwop,Buop; Czop,Dzwop,Dzuop; Cyop,Dywop,Dyuop};  % DJ (CC), 10/30/2026
is_cop = cellfun(@(P) isa(P,'copvar'),Pops);                                % DJ (CC), 10/30/2026
is_opv = cellfun(@(P) isa(P,'opvar') || isa(P,'opvar2d'),Pops);             % DJ (CC), 10/30/2026
if any(is_opv(:)) && ~any(is_cop(:))                                        % DJ (CC), 10/30/2026
    PIE = piess_opvar(Pops);                                                % DJ (CC), 10/30/2026
    return                                                                  % DJ (CC), 10/30/2026
end                                                                         % DJ (CC), 10/30/2026

% % % BEGIN construction from copvar operators, DJ (CC), 10/30/2026
% The construction from opvar operators that was here (DJ, 12/20/2024 -
% 10/29/2026) moved to piess_opvar, unchanged. The PIE is now built from
% the operators as they are, so that opvar and opvar2d operators among the
% copvar ones are converted by initialize, which also checks the spaces.
fnames = {'T','Tw','Tu'; 'A','Bw','Bu'; 'Cz','Dzw','Dzu'; 'Cy','Dyw','Dyu'};
% Type of component in the rows (x,x,z,y) and columns (x,w,u) of Pops.
rtype = repmat([1;1;2;3],1,3);
ctype = repmat(1:3,4,1);
is_mat = cellfun(@(P) (isa(P,'double') || isa(P,'polynomial')) && ~isempty(P),Pops);

% Registry of the variables of all copvar operators, and their domain.
cop = Pops(is_cop);
vars_all = cellfun(@(P) P.vars,cop(:)','UniformOutput',false);
dom_all = cellfun(@(P) P.dom,cop(:),'UniformOutput',false);
[vars,idx] = unique(reshape([cell(1,0),vars_all{:}],1,[]));
dom = cat(1,zeros(0,2),dom_all{:});
meta = struct('vars',{reshape(vars,1,[])},'dom',dom(idx,:));

% Convert matrices to copvar, between the spaces of the other operators.
spc = @(D,s) arrayfun(@(k) D.vars(D.(['space_',s])(k,:)),1:numel(D.(['dim_',s])),'UniformOutput',false);
for k = reshape(find(is_mat),1,[])
    Mat = double(Pops{k});
    dims = struct('out',size(Mat,1),'in',size(Mat,2));
    spaces = struct('out',{{{}}},'in',{{{}}});
    idr = find(is_cop & rtype==rtype(k),1);
    idc = find(is_cop & ctype==ctype(k),1);
    if ~isempty(idr)
        dims.out = Pops{idr}.dim_out;       spaces.out = spc(Pops{idr},'out');
    end
    if ~isempty(idc)
        dims.in = Pops{idc}.dim_in;         spaces.in = spc(Pops{idc},'in');
    end
    Pops{k} = mat2copvar_sop(Mat,dims,spaces,meta);
end

PIE = pie_struct();
for k = 1:numel(Pops)
    PIE.(fnames{k}) = Pops{k};
end
PIE = initialize(PIE);
% % % END construction from copvar operators, DJ (CC), 10/30/2026

end

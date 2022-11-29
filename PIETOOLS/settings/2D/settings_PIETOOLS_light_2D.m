%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings_PIETOOLS_light_2D.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function settings = settings_PIETOOLS_light_2D()
% Settings file for less challenging LPIs in two spatial variables.
%
% Complexity/accuracy of LPIs increases with the settings files as:
% extreme < stripped < light < heavy < veryheavy

%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP, SS, DJ - 02_21_2022

% Check if inputs are provided
if nargin>0
    error(['The settings functions_light does not take inputs.',...
        ' Please adjust the output settings struct, or adapt and run (a copy of) the "custom" settings file',...
        ' to enforce specific settings.'])
end
settings.is2D = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options for the Lyapunov Function V = <x,Px>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Set degrees of the monomials 
% % [Zx^o; Zx^a; Zx^b; 
% %  Zy^o; Zy^a; Zy^b; 
% %  Z2^oo; Z2^ao; Z2^bo; Z2^oa; Z2^ob; Z2^aa; Z2^ba; Z2^ab; Z2^bb];
% % defining the 2D PI operator P = Z'*Q*Z where Q>0 is positive definite.
% % See the "poslpivar_2d" header for more details.

% NOTE: The resulting operator is comprised of squares of these monomials
% --> degrees will doubled.
% NOTE: Unless the 2D PDE is coupled to 1D PDEs, degrees dx and dy will not
% contribute.

% Monomials in just the first spatial variable x
dx = {1;            % Degree of Zx^o(x) in x
      [1;1;2];      % Degree of Zx^a(x,tt) in x, tt and (x*tt)
      [1;1;2]};     % Degree of Zx^b(x,tt) in x, tt and (x*tt)
  
% Monomials in just the second spatial variable y
dy = {1,...         % Degree of Zy^o(y) in y
      [1,1,2],...   % Degree of Zy^a(y,nu) in y, nu and (y*nu)
      [1,1,2]};     % Degree of Zy^b(y,nu) in y, nu and (y*nu)
  
% Monomials in both spatial variables
d2 = {[0,1;1,2],...             % Degree of Z2^oo(x,y) in [-, y; x, x*y]
      [0,1,1,2;1,1,1,2],...     % Degree of Z2^oa(x,y,nu) in [-, y, nu, nu*y; x, x*y, x*nu, x*y*nu]
      [0,1,1,2;1,1,1,2];...     % Degree of Z2^ob(x,y,nu) in [-, y, nu, nu*y; x, x*y, x*nu, x*y*nu]
      [0,1;1,1;1,1;2,2],...     % Degree of Z2^ao(x,y,tt) in [-, y; x, x*y; tt, y*tt; x*tt, x*y*tt]
      [0,1,1,2;                 % Degree of Z2^aa(x,y,tt,nu) in [-,    y,      nu,      y*nu;
       1,2,2,2;                 %                                x,    x*y,    x*nu,    x*y*nu;
       1,2,2,2;                 %                                tt,   y*tt,   tt*nu,   y*tt*nu;
       2,2,2,2],...             %                                x*tt, x*y*tt, x*tt*nu, x*y*tt*nu]
      [0,1,1,2; 1,2,2,2; 1,2,2,2; 2,2,2,2];...      % Degree of Z2^ab(x,y,tt,nu)
      [0,1;1,1;1,1;2,2],...     % Degree of Z2^bo(x,y,tt) in [-, y; x, x*y; tt, y*tt; x*tt, x*y*tt]
      [0,1,1,2; 1,2,2,2; 1,2,2,2; 2,2,2,2],...      % Degree of Z2^ba(x,y,tt,nu)
      [0,1,1,2; 1,2,2,2; 1,2,2,2; 2,2,2,2]};        % Degree of Z2^bb(x,y,tt,nu)

settings.LF_deg.dx = dx;
settings.LF_deg.dy = dy;
settings.LF_deg.d2 = d2;
  
% % Should we exclude any parameters (monomials), or require separability?
settings.LF_opts.sep = [0;  % Set 1 to enforce Rxx{2}=Rxx{3}
                        0;  % Set 1 to enforce Ryy{2}=Ryy{3}
                        1;  % Set 1 to enforce R22{2,1}=R22{3,1}
                        1;  % Set 1 to enforce R22{1,2}=R22{1,2}
                        1;  % Set 1 to enforce R22{2,2}=R22{3,2}, R22{2,3}=R22{3,3}
                        1]; % Set 1 to enforce R22{2,2}=R22{2,3}, R22{3,2}=R22{3,3}
                    
% Set "exclude" value j equal to 1 to set Z(j)=0 in definition of P=Z'*Q*Z,
% where Z is defined by the 16 monomials
% [1;
%  Zx^o; Zx^a; Zx^b; 
%  Zy^o; Zy^a; Zy^b; 
%  Z2^oo; Z2^ao; Z2^bo; Z2^oa; Z2^ob; Z2^aa; Z2^ba; Z2^ab; Z2^bb];
settings.LF_opts.exclude = [0,...
                            0,0,0,...
                            0,0,0,...
                            0,0,0,0,0,1,1,1,1]; 

% % Should we add a Psatz term to require only local feasibility?
settings.LF_opts.psatz = 0;     % Do not enforce Psatz in initial operator P
settings.LF_use_psatz = 0;      % But perhaps add an operator P=P+P2 where P2 does use Psatz

% % We distinguish two psatz options:
% 1. Require P2>0 on square domain (x,y)\in[a,b]x[c,d]
settings.LF_opts_psatz{1}.psatz = 1;
% Otherwise, use the same settings as for P
settings.LF_opts_psatz{1}.sep = zeros(1,6);
settings.LF_opts_psatz{1}.exclude = zeros(1,16);
settings.LF_deg_psatz{1} = settings.LF_deg;

% 2. Require P2>0 on circle (x,y)\in {(x-0.5*(a+b))^2 + (y-0.5*(c+d))^2 <= r^2}
settings.LF_opts_psatz{2}.psatz = 2;
% Otherwise, use the same settings as for P
settings.LF_opts_psatz{2}.sep = zeros(1,6);
settings.LF_opts_psatz{2}.exclude = zeros(1,16);
settings.LF_deg_psatz{2} = settings.LF_deg;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options for the indeterminate PI variable Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In controller and estimator synthesis, the LPI depends on an additional,
% indeterminate variable Z, related to the optimal controller/estimator. We
% define this variable using the same monomials as for positive variables,
% but the indeterminate variable is not quadratic in these monomials.

% Monomials in just the first spatial variable x
Zop_dx = {2*dx{1};   2*dx{2};   2*dx{3}};
  
% Monomials in just the second spatial variable y
Zop_dy = {2*dy{1},   2*dy{2},   2*dy{3}};
  
% Monomials in both spatial variables
Zop_d2 = {2*d2{1,1},   2*d2{1,2},   2*d2{1,3};
          2*d2{2,1},   2*d2{2,2},   2*d2{2,3};
          2*d2{3,1},   2*d2{3,2},   2*d2{3,3}};
      
settings.Zop_deg.dx = Zop_dx;
settings.Zop_deg.dy = Zop_dy;
settings.Zop_deg.d2 = Zop_d2;

% For this operator, we can also enforce separability
settings.Zop_opts.sep = [0;  % Set 1 to enforce Rxx{2}=Rxx{3}
                         0;  % Set 1 to enforce Ryy{2}=Ryy{3}
                         0;  % Set 1 to enforce R22{2,1}=R22{3,1}
                         0;  % Set 1 to enforce R22{1,2}=R22{1,2}
                         0]; % Set 1 to enforce R22{2,2}=R22{3,2}=R22{2,3}=R22{3,3}
             
% We can require the parameters defining Zop to not vary in space, setting
% all degrees dx, dy and d2 to zero (independent of the values specified
% above).
settings.Zop_opts.isscalar = 0;     % Set to 1 for constant parameters

% We can also require Zop to be a multiplier operator, excluding any
% integrator terms. 
settings.Zop_opts.ismultiplier = 0; % Set to 1 to exclude integrators



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options for the LPI inequality constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We can enforce the inequality in two ways: using sosineq, or using soseq

% Establish if the user has already requested sosineq to be used
if evalin('base','exist(''sosineq_on'',''var'')') 
    use_sosineq = evalin('base','sosineq_on');
elseif evalin('base','exist(''use_sosineq'',''var'')') 
    use_sosineq = evalin('base','use_sosineq');
else
    use_sosineq = 1;
end
settings.use_sosineq = use_sosineq;

if use_sosineq
% Using sosineq, a constraint Q>=0 is enforced by building a new operator
% D>=0, and requiring Q==D. In doing so, the necessary degrees of the
% monomials defining D are established from the parameters defining Q.
% We have almost no freedom in constructing the operator D.    
    % We can, however, choose whether or not to add a Psatz term to D,
    settings.ineq_opts.psatz = 0;
    
else
% Useing soseq, a constraint Q==D for D>=0 is also enforced, but now we 
% define the operator D ourselves in the executive script. The degrees of D
% are determined by taking the degrees of Q, and adding some offset "Dup".

    % % Since we build D ourselves, we can set all defining options
    % Start with the degrees, increasing them slightly compared to those of Q
    
    % Monomials in just the first spatial variable x
    Dupx = 1;
    eq_dx = {Dupx+dx{1};   Dupx+dx{2};   Dupx+dx{3}};
    
    % Monomials in just the second spatial variable y
    Dupy = 1;
    eq_dy = {Dupy+dy{1},   Dupy+dy{2},   Dupy+dy{3}};
    
    % Monomials in both spatial variables
    Dup2 = 1;
    eq_d2 = {Dup2+d2{1,1},   Dup2+d2{1,2},   Dup2+d2{1,3};
             Dup2+d2{2,1},   Dup2+d2{2,2},   Dup2+d2{2,3};
             Dup2+d2{3,1},   Dup2+d2{3,2},   Dup2+d2{3,3}};
    
    settings.eq_deg.dx = eq_dx;
    settings.eq_deg.dy = eq_dy;
    settings.eq_deg.d2 = eq_d2;
         
    % Will we enforce any conditions on D?
    settings.eq_opts.psatz = 0;
    settings.eq_opts.sep = zeros(1,6);
    settings.eq_opts.exclude = zeros(1,16);    
    
    % Will we add a Psatz term?
    settings.eq_use_psatz = [0;2];          	% Set to 1 or 2 or [1;2] to use psatz
    for j=1:length(settings.eq_use_psatz)
        settings.eq_opts_psatz{j}.psatz = settings.eq_use_psatz(j);       
        settings.eq_opts_psatz{j}.sep = zeros(1,6);
        settings.eq_opts_psatz{j}.exclude = zeros(1,16);
        
        % Monomials in just the first spatial variable x
        eq_dx_psatz = {Dupx-1+dx{1};   Dupx-1+dx{2};   Dupx-1+dx{3}};
        eq_dy_psatz = {Dupy-1+dy{1},   Dupy-1+dy{2},   Dupy-1+dy{3}};
        eq_d2_psatz = {Dup2-1+d2{1,1},   Dup2-1+d2{1,2},   Dup2-1+d2{1,3};
                       Dup2-1+d2{2,1},   Dup2-1+d2{2,2},   Dup2-1+d2{2,3};
                       Dup2-1+d2{3,1},   Dup2-1+d2{3,2},   Dup2-1+d2{3,3}};

        settings.eq_deg_psatz{j}.dx = eq_dx_psatz;
        settings.eq_deg_psatz{j}.dy = eq_dy_psatz;
        settings.eq_deg_psatz{j}.d2 = eq_d2_psatz;
    end
end

end
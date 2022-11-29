function output = PIETOOLS_H2_norm_o(PIE, settings,varargin)
% This function solves the minimization problem to obtain the H2 norm of a linear distributed parameter
% system using the observability gramian approach and PIE framework. For
% the feasibility test an additional input with a assigned value to the
% norm is required.
% inputs: (mandatory)
%   PIE : PIE structure of the corresponding system
%   settings : options related to PIETOOLS
% output: structure with the following fields
%   h2= computed H2 norm for SDP problem or optional input value for feasibility tests;
%   W= Observability gramian, a positive PI operator;
%  prog= sum of squares program information;
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

   if PIE.dim==2
    % Call the 2D version of the executive.
      if nargin==2
            output = PIETOOLS_H2_norm_2D(PIE,settings);
      elseif nargin>2
           options=struct(varargin{:});              
            output = PIETOOLS_H2_norm_2D(PIE,settings,options);
      else
          disp('Incorrect number of inputs, try again \n')
          return
      end

      return

    end
    
    dd1 = settings.dd1;
    dd12 = settings.dd12;
    sos_opts = settings.sos_opts;
    options1 = settings.options1;
    options12 = settings.options12;
    override1 = settings.override1;
    eppos = settings.eppos;
    epneg = settings.epneg;
    eppos2 = settings.eppos2;
    ddZ = settings.ddZ;
    sosineq_on = settings.sosineq_on;
    if sosineq_on
        opts = settings.opts;
    else
        override2 = settings.override2;
        options2 = settings.options2;
        options3 = settings.options3;
        dd2 = settings.dd2;
        dd3 = settings.dd3;
    end
    
    % Dumping relevant 4-PI operators to the workspace 
    Aop=PIE.A;
    Top=PIE.T;
    B1op=PIE.B1;    %TB1op = PIE.Tw;
    C1op=PIE.C1;
    %D11op=PIE.D11;
    
    fprintf('\n --- Searching for H2 norm bound using the observability gramian --- \n')
    % Declare an SOS program and initialize domain and opvar spaces
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    varlist = [Aop.var1; Aop.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
    prog = sosprogram(varlist);      % Initialize the program structure
    X=Aop.I;                         % retrieve the domain from Aop
    nx1=Aop.dim(1,1);                % retrieve the number of ODE states from Aop
    nx2=Aop.dim(2,1);                % retrieve the number of distributed states from Aop
    nw=B1op.dim(1,2);                % retrieve the number of real-valued disturbances
    nz=C1op.dim(1,1);                % retrieve the number of real-valued regulated outputs
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % if the user wants to calculate the norm, define objective function
    % otherwise, the norm is input to the function.
    if nargin<=2 || ~isfield(options,'h2')
        dpvar gam;
        prog = sosdecvar(prog, gam); %this sets gam = gamma as decision var
        prog = sosineq(prog, gam); %this ensures gamma is lower bounded
        prog = sossetobj(prog, gam); %this minimizes gamma, comment for feasibility test
    else
        gam = options.h2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1: declare the posopvar variable, Pop, which defines the storage 
    % function candidate
    disp('- Declaring Gramian using specified options...');
    
    [prog, P1op] = poslpivar(prog, [nx1 ,nx2],X,dd1,options1);
    
    if override1~=1
        [prog, P2op] = poslpivar(prog, [nx1 ,nx2],X,dd12,options12);
        Wop=P1op+P2op;
    else
        Wop=P1op;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 2: Using the observability gramian
    
    disp('- Constructing the Negativity Constraint...');
    
    Dop =  (Aop'*Wop)*Top+Top'*(Wop*Aop)+C1op'*C1op;
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 3: Impose Negativity Constraint. There are two methods, depending on
    % the options chosen
    %
    disp('- Enforcing the Inequalities Constraints...');
    if sosineq_on
        disp('  - Using lpi_ineq...');
        prog = lpi_ineq(prog,-Dop,opts);
    else
        disp('  - Using an Equality constraint...');
        [prog, De1op] = poslpivar(prog, [nx1, nx2],X,dd2,options2);
        
        if override2~=1
            [prog, De2op] = poslpivar(prog,[nx1, nx2],X, dd3,options3);
            Deop=De1op+De2op;
        else
            Deop=De1op;
        end
        prog = lpi_eq(prog,Deop+Dop); %Dop=-Deop
    end
    
    tempObj = B1op'*Wop*B1op;
    tempMat = tempObj.P;
    traceVal=0;
    for idx = 1:size(tempMat,1)
        traceVal = traceVal+tempMat(idx,idx);
    end
    % ensuring scalar inequality gam>trace
    prog = sosineq(prog, gam-traceVal);
    
    %solving the sos program
    disp('- Solving the LPI using the specified SDP solver...');
    prog = sossolve(prog,sos_opts); 
    
    % disp('The H2 norm of the given system is upper bounded by:')
    % if ~isreal(gam)
    %     disp(double(sosgetsol(prog,gam))); % check the Hinf norm, if the solved successfully
    % else 
    %     disp(gam);
    % end
     P = getsol_lpivar(prog,Wop);
    % gam = double(sosgetsol(prog,gam));
    % end
   if nargin<=2 || ~isfield(options,'h2')
            gam = sqrt(double(sosgetsol(prog,gam)));
            disp('The H2 norm of the given system is upper bounded by:')
             disp(gam);
   end
%      gam = sqrt(double(sosgetsol(prog,gam)));
%             disp('The H2 norm of the given system is upper bounded by:')
%              disp(gam);
    %% set outputs
    output.h2=gam;
    output.W=P;
    output.prog=prog;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

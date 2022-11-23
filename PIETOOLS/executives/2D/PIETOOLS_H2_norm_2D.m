function output = PIETOOLS_H2_norm_2D(PIE, settings,options)
% STEP 0: Extract LPI settings and necessary PI operators
    

%     if nargin==1
%         settings = settings_PIETOOLS_light_2D;
%         settings.sos_opts.simplify = 1;         % Use psimplify
%         settings.eppos = [1e-4; 1e-6; 1e-6; 1e-6];    % Positivity of Lyapunov Function
%         settings.epneg = 0*1e-5;                % Negativity of derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
%     elseif ~isfield(settings,'is2D') || ~settings.is2D
        % Extract 2D settings.
        sos_opts = settings.sos_opts;
        settings = settings.settings_2d;
        settings.sos_opts = sos_opts;
   % end
    
    % Extract options for sossolve
    if ~isfield(settings,'sos_opts')
        sos_opts.simplify = 1;         % Use psimplify
        sos_opts.solver = 'sedumi';
        % % Other optional SDP solvers supported by PIETOOLS
        % settings.sos_opts.solver = 'mosek';
        % settings.sos_opts.solver = 'sdpnalplus';
        % settings.sos_opts.solver = 'sdpt3';
    else
       sos_opts = settings.sos_opts; 
    end
    
    % Extract LF positivity and negativity conditions
    if ~isfield(settings,'eppos')
        eppos = [1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
                  1e-6;1e-6;
                  1e-6];    % Positivity of Lyapunov Function with respect to 2D spatially distributed states
    else
        eppos = settings.eppos;
    end
    if ~isfield(settings,'epneg')
        epneg = 0;         % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
    else
        epneg = settings.epneg;
    end
    
    % Extract settings defining operator P parameterizing LF V=<x,Px>
    LF_deg = settings.LF_deg;
    LF_opts = settings.LF_opts;
    
    % Does P include a psatz term? If so, what degrees and options?
    LF_use_psatz = settings.LF_use_psatz;
    LF_deg_psatz = extract_psatz_deg(settings.LF_deg_psatz,LF_use_psatz);
    LF_opts_psatz = extract_psatz_opts(settings.LF_opts_psatz,LF_use_psatz);
    
    % Do we use lpi_ineq or not? Extract the appropriate options
    use_lpi_ineq = settings.use_sosineq;
    if use_lpi_ineq
        ineq_opts = settings.ineq_opts;
    else
        eq_opts = settings.eq_opts;
        eq_deg = settings.eq_deg;
        
        eq_use_psatz = settings.eq_use_psatz;
        eq_deg_psatz = extract_psatz_deg(settings.eq_deg_psatz,eq_use_psatz);
        eq_opts_psatz = extract_psatz_opts(settings.eq_opts_psatz,eq_use_psatz);
    end
    
    % Extract the relevant PI operators
    Aop=PIE.A;
    Top=PIE.T;
    B1op=PIE.B1;    %TB1op = PIE.Tw;
    C1op=PIE.C1;
    %D11op=PIE.D11;
    % Other relevant information
    vars = PIE.vars;                % Spatial variables of the PIE
    dom = PIE.dom;                  % Spatial domain of the PIE
    np_op = Aop.dim(:,1);           % Dimensions RxL2[x]xL2[y]xL2[x,y] of the PDE state
    
    fprintf('\n --- Searching for H2 norm bound using the observability gramian --- \n')
    % Declare an SOS program and initialize domain and opvar spaces
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    prog = sosprogram(vars(:));      % Initialize the program structure
    prog.vartable = [vars(:,1).varname; vars(:,2).varname]; % Make sure the variables are in the right order.
    
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
    %
    % Alternatively, the above 3 commands may be commented and a specific gain
    % test specified by defining a specific desired value of gamma. This
    % results in a feasibility test instead of an optimization problem.
    % gamma = 1000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 1: declare the posopvar variable, Pop, which defines the storage 
    % function candidate
    disp('- Declaring Gramian using specified options...');
    
    % Initialize an operator which is positive semidefinite everywhere
    [prog, P1op] = poslpivar_2d(prog,np_op,dom,LF_deg,LF_opts);
    %[prog, P1op] = poslpivar(prog, [nx1 ,nx2],X,dd1,options1);
    
    % Add an additional term with psatz multiplier if requested
    for j=1:length(LF_use_psatz)
        if LF_use_psatz(j)~=0
            [prog, P2op] = poslpivar_2d(prog,np_op,dom,LF_deg_psatz{j},LF_opts_psatz{j});
            Wop = P1op + P2op;
        end
    end
    
    % Ensure strict positivity of the operator
    if ~all(eppos==0)
        Ip = blkdiag(eppos(1)*eye(np_op(1)),zeros(np_op(2)),zeros(np_op(3)),eppos(4)*eye(np_op(4)));
        Iop = opvar2d(Ip,[np_op,np_op],dom,vars);
        Wop = P1op + Iop;
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
    disp('- Enforcing the Inequality Constraints...');

     if use_lpi_ineq
        disp('   - Using lpi_ineq...');
        prog = lpi_ineq_2d(prog,-Dop,ineq_opts);
    else
        disp('   - Using an equality constraint...');
        [prog, Deop]  =poslpivar_2d(prog,np_op,dom,eq_deg,eq_opts);
        
        % Introduce the psatz term.
        for j=1:length(eq_use_psatz)
            if eq_use_psatz(j)~=0
                eq_opts_psatz{j}.exclude(8:12) = eq_opts_psatz_exclude(8:12);
                [prog, De2op] = poslpivar_2d(prog,np_op+nw_op+nz_op,dom, eq_deg_psatz{j},eq_opts_psatz{j});
                Deop = Deop+De2op;
            end
    
        end

        prog = lpi_eq_2d(prog,Deop+Dop); %Dop=-Deop
     end
    
    tempObj = B1op'*Wop*B1op;
    tempMat = tempObj.R00;
    traceVal=0;
    for idx = 1:size(tempMat,1)
        traceVal = traceVal+tempMat(idx,idx);
    end
    % traceVal>=gam
    prog = sosineq(prog, gam-traceVal);
    
    
    %solving the sos program
    disp('- Solving the LPI using the specified SDP solver...');
    prog = sossolve(prog,sos_opts); 
    
     P = getsol_lpivar(prog,Wop);
    % gam = double(sosgetsol(prog,gam));
    % end
     if nargin<=2 || ~isfield(options,'h2')
            gam = sqrt(double(sosgetsol(prog,gam)));
            disp('The H2 norm of the given system is upper bounded by:')
             disp(gam);
    end
    %% set outputs
    output.h2=gam;
    output.W=P;
    output.prog=prog;

end

function outcell = extract_psatz_deg(incell,use_psatz)
% Check the number of elements of incell against those of use_psatz to
% determine if a cell element has been defined for each psatz term.

    if all(use_psatz==0)
        outcell = {};
        return
    end
    if isa(incell,'struct')
        for j=1:length(use_psatz)
            outcell{j} = incell;
        end
    elseif numel(incell)==1
        for j=1:length(use_psatz)
            outcell{j} = incell{1};
        end
    elseif numel(incell)==length(use_psatz)
        outcell = incell;
    elseif numel(incell)>=max(use_psatz)
        outcell = incell(use_psatz);
    else
        error('For each element of ''use_psatz'', a ''deg'' field should be defined.')
    end

end
function outcell = extract_psatz_opts(incell,use_psatz)
% Check the number of elements of incell against those of use_psatz to
% determine if a cell element has been defined for each psatz term.
    
    if all(use_psatz==0)
        outcell = {};
        return
    end
    if isa(incell,'struct')
        for j=1:length(use_psatz)
            outcell{j} = incell;
        end
    elseif numel(incell)==1
        for j=1:length(use_psatz)
            outcell{j} = incell{1};
        end
    elseif numel(incell)==length(use_psatz)
        outcell = incell;
    elseif numel(incell)>=max(use_psatz)
        outcell = incell(use_psatz);
    else
        error('For each element of ''use_psatz'', an ''opts'' field should be defined.')
    end

end
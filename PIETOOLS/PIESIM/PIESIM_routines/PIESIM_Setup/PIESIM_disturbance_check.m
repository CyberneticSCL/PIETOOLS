%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_disturbance_check.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uinput = PIESIM_disturbance_check(varargin)
% 1) Check if disturbances are provided in a correct format 
% 2) Split diturbances into a spatial and temporal parts for simulation
% Inputs:
% 1) varargin{1}: structure - PDE or PIE 
% 2) varargin{2}: uinput - user-defined boundary inputs, forcing and initial
% conditions. If empty or incomplete, will be set to default values 
% 3) varargin{3}: psize - size of the problem
% All inputs are REQUIRED
% Outputs: uinput -user-defined boundary inputs, forcing and initial
% conditions. If empty or incomplete, will be set to default values 

% Initial coding: YP 1/4/2026

syms st sx sy;

nargin=length(varargin);

structure=varargin{1};
uinput=varargin{2};
psize=varargin{3};

if isa(structure,'pie_struct') || isfield(structure,'T')
    type='PIE';
    PIE=structure;
elseif isa(structure,'pde_struct') || isa(structure,'sys')||isfield(structure,'BC') || isfield(structure,'X')
    type='PDE';
    PDE=structure;
else
    error('Data structure must be provided to the disturbance check');
end

 % Disturbance check for PDEs

nw=psize.nw;
nw0=psize.nw0;
nu=psize.nu;
nu0=psize.nu0;

  if (nu>0 & ~isfield(uinput,'u') )
            if ~isfield(uinput,'control')
                disp('Signal type for control inputs is not defined. Defaulting to zero');
                       uinput.u{1:nu}=0;
            else
            disp('Control inputs are not defined. Defaulting to a user-specified signal type');
                switch uinput.control
                    case 'constant'
                        [uinput.u{1:nu}]=deal(1+0*st);
                    case 'sin'
                        [uinput.u{1:nu}]=deal(sin(st));
                    case 'sinc'
                        [uinput.u{1:nu}]=deal(sinc(st));
                    otherwise
                        disp('The specified signal type for control inputs is not supported. Defaulting to zero');
                       uinput.u{1:nu}=0;
                end
            end
    end


    % From here we assume that uinput.u exists
         if (nu>0 | isfield(uinput,'u'))

         % Make sure uinput.u is a cell array
            if ~iscell(uinput.u)
            uinput.u=num2cell(uinput.u);
            end

            % Check control inputs size and type
            
                 if (size(uinput.u,2)<nu)
                    disp('Warning: Number of provided u inputs is less than nu.');
                    disp('Defaulting the rest of u inputs to zero');
                    uinput.u(size(uinput.u,2)+1:nu)={0};
                elseif (size(uinput.u,2)>nu)
                    disp('Warning: Number of provided u inputs is  greater than nu.');
                    disp('The rest of u inputs will be ignored');
                    uinput.u=uinput.u(1:nu);
                 end
         end % nu>0

         for k=1:nu
     % Check that the input is of appropriate type 
        if ~isdouble(uinput.u{k}) & ~isa(uinput.u{k}, 'sym') 
            fprintf('Warning: uinput.u for the input %d is neither ''double'' nor ''sym''. Defaulting to zero.\n', k);
            uinput.u{k}=0;
        end
        % Check size for double
         if isdouble(uinput.u{k})

           if min(size(uinput.u{k}))==1 & max(size(uinput.u{k}))>1 
            fprintf('Warning: uinput.u for the input %d needs to be or size 2 x nstep. Defaulting to zero. \n', k)
             uinput.u{k}=0;
             % Make sure control inputs are of size 2 x nstep
           elseif size(uinput.u{k},1)>size(uinput.u{k},2)
               uinput.u{k}=transpose(uinput.u{k});
           end

         end % ~isdouble(uinput.u{k})
        end % k

        % Check disturbance setup

    if (nw>0 & ~isfield(uinput,'w') )
            if ~isfield(uinput,'dist')
                disp('Signal type for disturbances is not defined. Defaulting to zero');
                       [uinput.w{1:nw}] = deal(0);
            else
            disp('Disturbances are not defined. Defaulting to a user-specified signal type');
                switch uinput.dist
                    case 'constant'
                        [uinput.w{1:nw}]=deal(1+0*st);
                    case 'sin'
                        [uinput.w{1:nw}]=deal(sin(st));
                    case 'sinc'
                        [uinput.w{1:nw}]=deal(sinc(st));
                    otherwise
                        disp('The specified signal type for disturbances is not supported. Defaulting to zero');
                       [uinput.w{1:nw}] = deal(0);
                end
            end
    end

     % From here we assume that uinput.w exists
         if (nw>0 | isfield(uinput,'w'))

              % Make sure uinput.w is a cell array
            if ~iscell(uinput.w)
            uinput.w=num2cell(uinput.w);
            end

                if (size(uinput.w,2)<nw)
                    disp('Warning: Number of provided w inputs is less than nw');
                    disp('Defaulting the rest of w inputs to zero');
                    uinput.w(size(uinput.w,2)+1:nw)={0};
                elseif (size(uinput.w,2)>nw)
                    disp('Warning: Number of provided w inputs is greater than nw');
                    disp('Extra disturbances will be ignored');
                    uinput.w=uinput.w(1:nw);
                end

         end % nw>0

    for k=1:nw
     % Check that the disturbance is of appropriate type 
        if ~isdouble(uinput.w{k}) & ~isa(uinput.w{k}, 'sym')
            fprintf('Warning: uinput.w for the input %d is neither ''double'' nor ''sym''. Defaulting to zero.\n', k);
            uinput.w{k}=0;
        end
        % Check size for double
         if isdouble(uinput.w{k})

           if min(size(uinput.w{k}))==1 & max(size(uinput.w{k}))>1 
             fprintf('Warning: uinput.w for the input %d needs to be of size 2 x nstep. Defaulting to zero. \n', k)
             uinput.w{k}=0;
             % Make sure disturbances are of size 2 x nstep
           elseif size(uinput.w{k},1)>size(uinput.w{k},2)
               uinput.w{k}=transpose(uinput.w{k});
           end

         end % isdouble(uinput.w{k})
        end % k
%   ------------------------------------------------------------------------------------------------
%   Disturbance and control inputs are now checked. 
%   Proceed to reordering and splitting.
%   ------------------------------------------------------------------------------------------------

         if type=='PDE' 
        [PDE,w_order] = reorder_comps(PDE,'w'); % Reorder disturbances in increasing order of differentiability
        [PDE,u_order] = reorder_comps(PDE,'u'); % Reorder control inputs in increasing order of differentiability 
             
         % Append this to uinput for future use
         uinput.w_order=w_order;
         uinput.u_order=u_order;

         if~issorted(u_order)
         uinput.u = uinput.u(u_order); % Reorder control inputs if needed
         end

         if~issorted(w_order)
         uinput.w = uinput.w(w_order); % Reorder disturbances if needed
         end

         end   % PDE

    % Split disturbances into spatial and temporal parts
    for k=1:nw
        uinput.wsep{k}=true;
         if (k<=nw0)
                dom=[];
         else
            if (type=='PDE')
        dom=PDE.w{k,1}.dom;
            else
        dom=PIE.dom;
            end
         end % if k<=nw0
    if(~isdouble(uinput.w{k}))
  entry=uinput.w{k};
    if(has(entry,sx)|has(entry,sy))
        split = cell2sym(children(entry));
        temp=split(has(split,st));
        temp=prod(temp);

        if (isempty(temp))
        uinput.wspace{k}=uinput.w{k};
        uinput.w{k}=1+0*st;
        elseif(has(temp,sx)|has(temp,sy)) % temp non-empty but split unsuccessful
                 uinput.wsep{k}=false;
                 uinput.wspace{k}=sym(0);
                 entry=0;
        else % temp non-empty but split successful
         uinput.wspace{k}=uinput.w{k}/temp;
         uinput.w{k}=temp;
        end % if (isempty(temp))
        
        % Check if disturbance is in a separable form
        if has(uinput.wspace{k},st)
                  uinput.wsep{k}=false;
                  entry=0;
        end % has(uinput.wspace{k},st)

          % Check if disturbance is specified correctly
        if (isempty(dom))
            disp('Warning: spatial contribution of disturbance is provided for a finite-dimensional disturbance. Spatial part will be ignored.');
            if (isdouble(entry))
            uinput.w{k}=entry*temp;
            else
            uinput.w{k}=coeffs(entry)*temp;
            end
            uinput.wspace{k}=sym(0);
        end
       
    else % entry does not have sx or sy - disturbance is symbolic but space-independent
      if (isempty(dom))
        uinput.wspace{k}=sym(0);
        else
        uinput.wspace{k}=sym(1);
        end
    end %  (has(entry,sx)|has(entry,sy)) 

    else % uinput.w{k} is double
        % Double
        if (isempty(dom))
        uinput.wspace{k}=sym(0);
        else
        uinput.wspace{k}=sym(1);  
        end
    end % (~isdouble(uinput.w{k}))

    end % for k

        % Split control inputs into spatial and temporal parts
    for k=1:nu
        uinput.usep{k}=true;
         if (k<=nu0)
                dom=[];
         else
            if (type=='PDE')
        dom=PDE.u{k,1}.dom;
            else
        dom=PIE.dom;
            end
         end % if k<=nu0
    if(~isdouble(uinput.u{k}))
  entry=uinput.u{k};
    if(has(entry,sx)|has(entry,sy))
        split = cell2sym(children(entry));
        temp=split(has(split,st));
        temp=prod(temp);

        if (isempty(temp))
        uinput.uspace{k}=uinput.u{k};
        uinput.u{k}=1+0*st;
        elseif(has(temp,sx)|has(temp,sy)) % temp non-empty but split unsuccessful
                 uinput.usep{k}=false;
                 uinput.uspace{k}=sym(0);
                 entry=0;
        else % temp non-empty but split successful
         uinput.uspace{k}=uinput.u{k}/temp;
         uinput.u{k}=temp;
        end % if (isempty(temp))
        
        % Check if disturbance is in a separable form
        if has(uinput.uspace{k},st)
                  uinput.usep{k}=false;
                  entry=0;
        end % has(uinput.wspace{k},st)

          % Check if disturbance is specified correctly
        if (isempty(dom))
            disp('Warning: spatial contribution of control input is provided for a finite-dimensional disturbance. Spatial part will be ignored.');
            if (isdouble(entry))
            uinput.u{k}=entry*temp;
            else
            uinput.u{k}=coeffs(entry)*temp;
            end
            uinput.uspace{k}=sym(0);
        end
       
    else % entry does not have sx or sy - disturbance is symbolic but space-independent
      if (isempty(dom))
        uinput.uspace{k}=sym(0);
        else
        uinput.uspace{k}=sym(1);
        end
    end %  (has(entry,sx)|has(entry,sy)) 

    else % uinput.w{k} is double
        % Double
        if (isempty(dom))
        uinput.uspace{k}=sym(0);
        else
        uinput.uspace{k}=sym(1);  
        end
    end % (~isdouble(uinput.u{k}))

    end % for k

 

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_DDF.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DDF_out=initialize_PIETOOLS_DDF(DDF)
% This a initialize/preprocessing routine for a DDF with K = numel(tau)
% delays,
%   \dot x(t) = A0 x(t) + B1 w(t) + B2 u(t) + Bv v(t)
%        z(t) = C1 x(t) + D11 w(t) + D12 u(t) + D1v v(t)
%        y(t) = C2 x(t) + D21 w(t) + D22 u(t) + D2v v(t)
%      r_i(t) = Cr{i} x(t) + Br1{i} w(t) + Br2{i} u(t) + Drv{i} v(t)
%        v(t) = sum_i Cv{i} r_i(t-tau(i)) + sum_i int_{-tau(i)}^0 Cvd{i}(s) r_i(t+s) ds
% The preprocessing stage checks which matrices have been defined and
% reads the dimensions nx, nz, ny, nr{i}, nv, nw, nu off them; undefined
% matrices are set to zero, defined ones are checked against these
% dimensions.
% undefined signals are set to length 0 with 0 value
%
% input: DDF - a partially-defined DDF data structure: 'struct' with the
%        field tau (1 x K 'double'; required), the 'double' matrices A0,
%        B1, B2, Bv, C1, D11, D12, D1v, C2, D21, D22, D2v, cells of at
%        most K 'double' matrices Cr, Br1, Br2, Drv, Cv and of 'double'
%        or 'polynomial' kernels Cvd (empty or missing elements are zero),
%        and optional fields nx, nz, ny, nw, nu that fix the dimensions.
%
% output: DDF_out - a fully-defined DDF data structure: 'struct' with the
%        fields tau, A0, B1, B2, C1, C2, D11, D12, D21, D22, Bv, D1v, D2v
%        and the 1 x K cells Cr, Br1, Br2, Drv, Cv, Cvd.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding MP - 10_01_2020
% MMP, 5_30_2021; added new data structure and made script a function
% MMP, 11/21/2026: Changes since PIETOOLS 2025: nx can be read off B1; a
%   coefficient with either dimension wrong throws the size error; help and
%   comments revised.


%%%%%%%%%%%%%%%%%%%%%%%%%
% Dumping the data structure to the workspace. Stop judging us!!!
names=fieldnames(DDF);
for i=1:numel(names)
    eval([names{i},'=DDF.',names{i},';'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking for undefined matrices
%%%%%%%%%%%%%%%%%%%%%%% Group 1 %%%%%%%%%%%%%%%%%%%%%%%

if ~exist('tau')
    error('you forgot to define the delays, tau')
end
nK=length(tau);    % number of delays

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Before Proceeding, we should check that the user has defined enough 
%%%% delays. For this we check the number of elements in each of the cell
%%%% arrays to make sure they are at most nK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('Cr','var')&&(numel(Cr)>nK)
    error('You have defined an element of Cr for which there is no matching value of delay.')
end

if exist('Br1','var')&&(numel(Br1)>nK)
    error('You have defined an element of Br1 for which there is no matching value of delay.')
end

if exist('Br2','var')&&(numel(Br2)>nK)
    error('You have defined an element of Br2 for which there is no matching value of delay.')
end

if exist('Drv','var')&&(numel(Drv)>nK)
    error('You have defined an element of Drv for which there is no matching value of delay.')
end

if exist('Cv','var')&&(numel(Cv)>nK)
    error('You have defined an element of Cv for which there is no matching value of delay.')
end

if exist('Cvd','var')&&(numel(Cvd)>nK)
    error('You have defined an element of Cvd for which there is no matching value of delay.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First determine all signal sizes
% nx, nz, ny, nri, nv, nw, nu
%%%%%%%%%%%%%%%%%%%%%%% Group 1 %%%%%%%%%%%%%%%%%%%%%%%
nstates=0;noregs=0;nomeas=0;nori=0;nov=0;now=0;nou=0; % toggles to determine if no states, inputs or outputs are present

% scanning for the number of states
if ~exist('nx','var')
    if exist('A0','var')||exist('B1','var')||exist('B2','var')||exist('Bv','var')
        % nx: rows of A0, B1, B2 or Bv, the first one defined
        if exist('A0','var')
            nx=size(A0,1);
        elseif exist('B1','var')
%           nx=size(B1,1,'var');                                            % MMP, 11/21/2026 (was)
            nx=size(B1,1);                                                  % MMP, 11/21/2026
        elseif exist('B2','var')
            nx=size(B2,1);
        elseif exist('Bv','var')
            nx=size(Bv,1);
        end
    else
        nx=0;
        warning('No state dynamics detected. There are no states')
        nostates=1;
    end
end

% With no states (nx = 0), C1, C2 and Cr{i} have no columns.

% scanning for the number of regulated outputs
if ~exist('nz','var')
    if exist('C1','var')||exist('D11','var')||exist('D12','var')||exist('D1v','var')
        % nz: rows of C1, D11, D12 or D1v, the first one defined
        if exist('C1','var')
            nz=size(C1,1);
        elseif exist('D11','var')
            nz=size(D11,1);
        elseif exist('D12','var')
            nz=size(D12,1);
        elseif exist('D1v','var')
            nz=size(D1v,1);
        end
    else
        disp('No regulated output channel detected')
        nz=0;
        noregs=1;
    end
end

% scanning for the number of sensed outputs
if ~exist('ny','var')
    if exist('C2','var')||exist('D21','var')||exist('D22','var')||exist('D2v','var')
        % ny: rows of C2, D21, D22 or D2v, the first one defined
        if exist('C2','var')
            ny=size(C2,1);
        elseif exist('D21','var')
            ny=size(D21,1);
        elseif exist('D22','var')
            ny=size(D22,1);
        elseif exist('D2v','var')
            ny=size(D2v,1);
        end
    else
        disp('No measured output channel detected')
        ny=0;
        nomeas=1;
    end
end

% Detect if there are unused delay channels
if exist('Cv','var')
    lcv=numel(Cv);
else
    lcv=0;
end
if exist('Cvd','var')
    lcvd=numel(Cvd);
else
    lcvd=0;
end

% No Cv and no Cvd defined: warning
if lcvd==0&&lcv==0
        warning(['Problem: The delayed channels have no outputs'])
        nv=0;
end    
nv=0;
% nv: rows of Cv{i}, else of Cvd{i}, of the last channel with an output
for i=1:nK
    if lcv>=i && ~isempty(Cv{i})
        nv=size(Cv{i},1);
    elseif lcvd>=i && ~isempty(Cvd{i})
        nv=size(Cvd{i},1);
    else
        warning(['Problem: Delayed channel ',int2str(i),' has no outputs'])
        nov=1;
    end
    
end



% Now we scan for the number of states in each delayed channel, r_i
% The inputs to channel r_i are contained in Cr{i}, Br1{i}, Br2{i}, Drv{i}
if exist('Cr','var')
    lcr=numel(Cr);
else
    lcr=0;
end
if exist('Br1','var')
    lbr1=numel(Br1);
else
    lbr1=0;
end
% Numbers of elements of Br2 and Drv
if exist('Br2','var')
    lbr2=numel(Br2);
else
    lbr2=0;
end
if exist('Drv','var')
    ldrv=numel(Drv);
else
    ldrv=0;
end

% nr{i}: rows of the first defined of Cr{i}, Br1{i}, Br2{i}, Drv{i}
for i=1:nK
    if lcr>=i&&~isempty(Cr{i})
        nr{i}=size(Cr{i},1);   % number of states
    elseif lbr1>=i&&~isempty(Br1{i})
        nr{i}=size(Br1{i},1);   % number of states
    elseif lbr2>=i&&~isempty(Br2{i})
        nr{i}=size(Br2{i},1);   % number of states
    elseif ldrv>=i&&~isempty(Drv{i})
        nr{i}=size(Drv{i},1);   % number of states
    else
        warning(['Problem: Delayed channel ',int2str(i),' has no inputs'])
        nr{i}=0;
    end
end





% scanning for the number of disturbance inputs
if ~exist('nw','var')
    if exist('B1','var')||exist('D11','var')||exist('D21','var')||exist('Br1','var')
        % nw: columns of B1, D11 or D21, the first one defined, else of the last element of Br1
        if exist('B1','var')
            nw=size(B1,2);
        elseif exist('D11','var')
            nw=size(D11,2);
        elseif exist('D21','var')
            nw=size(D21,2);
        elseif exist('Br1','var')
            if isempty(Br1{numel(Br1)})
                error('The last entry in Br1 is empty')
            end
            nw=size(Br1{numel(Br1)},2);
        end
    else
        % No disturbance coefficient: nw = 0
        disp('No disturbance channels detected')
        now=1;
        nw=0;
    end
end

% scanning for the number of control inputs
if ~exist('nu','var')
    if exist('B2','var')||exist('D12','var')||exist('D22','var')||exist('Br2','var')
        % nu: columns of B2, D12 or D22, the first one defined, else of the last element of Br2
        if exist('B2','var')
            nu=size(B2,2);
        elseif exist('D12','var')
            nu=size(D12,2);
        elseif exist('D22','var')
            nu=size(D22,2);
        elseif exist('Br2','var')
            if isempty(Br2{numel(Br2)})
                error('The last entry in Br2 is empty')
            end
            
            nu=size(Br2{numel(Br2)},2);
        end
    else
        % No control-input coefficient: nu = 0
        disp('No control input channels detected')
        nu=0;
        nou=1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Now that we know the number of inputs and outputs, we can check each
%%%% of the user-defined matrices. If they don't exist they are defaulted
%%%% to zero. If they exist, we verify they are of the correct dimension.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Group 1 %%%%%%%%%%%%%%%%%%%%%%
% A0: zero if undefined, else size checked (nx x nx)
if ~exist('A0','var')
    A0=zeros(nx,nx);
    disp('A0 undefined, defaulting to 0')
elseif (size(A0,1)~=nx)||(size(A0,2)~=nx)
    error('The dimensions of A0 are inconsistent with the detected number of inputs or states')
end


% B1: zero if undefined, else size checked (nx x nw)
if ~exist('B1','var')
    B1=zeros(nx,nw);
    disp('B1 undefined, defaulting to 0')
elseif (size(B1,1)~=nx)||(size(B1,2)~=nw)
    error('The dimensions of B1 are inconsistent with the detected number of inputs or states')
end

% B2: zero if undefined, else size checked (nx x nu)
if ~exist('B2','var')
    B2=zeros(nx,nu);
    disp('B2 undefined, defaulting to 0')
elseif (size(B2,1)~=nx)||(size(B2,2)~=nu)
    error('The dimensions of B2 are inconsistent with the detected number of inputs or states')
end

% Bv: zero if undefined, else size checked (nx x nv)
if ~exist('Bv','var')
    Bv=zeros(nx,nv);
    disp('Bv undefined, defaulting to 0')
elseif (size(Bv,1)~=nx)||(size(Bv,2)~=nv)
    error('The dimensions of Bv are inconsistent with the detected number of inputs or states')
end



%%%%%%%%%%%%%%%%%%%%%%% Group 2 %%%%%%%%%%%%%%%%%%%%%%%

for i=1:nK
    % Cr{i}, Br1{i}: zero where undefined or empty, else size checked
    % (nr{i} x nx, nr{i} x nw) in both dimensions                           % MMP, 11/21/2026
    if lcr<i||isempty(Cr{i})
        Cr{i}=zeros(nr{i},nx);
        disp(['Cr{',int2str(i),'} undefined, defaulting to 0'])
%   elseif size(Cr{i})~=[nr{i},nx]                                          % MMP, 11/21/2026 (was)
    elseif ~isequal(size(Cr{i}),[nr{i},nx])                                 % MMP, 11/21/2026
        error(['error in size of Cr',int2str(i)])
    end
    
    if lbr1<i||isempty(Br1{i})
        Br1{i}=zeros(nr{i},nw);
        disp(['Br1{',int2str(i),'} undefined, defaulting to 0'])
%   elseif size(Br1{i})~=[nr{i},nw]                                         % MMP, 11/21/2026 (was)
    elseif ~isequal(size(Br1{i}),[nr{i},nw])                                % MMP, 11/21/2026
        error(['error in size of Br1',int2str(i)])
    end
    
    % Br2{i}, Drv{i}: likewise (nr{i} x nu, nr{i} x nv)
    if lbr2<i||isempty(Br2{i})
        Br2{i}=zeros(nr{i},nu);
        disp(['Br2{',int2str(i),'} undefined, defaulting to 0'])
%   elseif size(Br2{i})~=[nr{i},nu]                                         % MMP, 11/21/2026 (was)
    elseif ~isequal(size(Br2{i}),[nr{i},nu])                                % MMP, 11/21/2026
        error(['error in size of Br2',int2str(i)])
    end
    
    if ldrv<i||isempty(Drv{i})
        Drv{i}=zeros(nr{i},nv);
        disp(['Drv{',int2str(i),'} undefined, defaulting to 0'])
%   elseif size(Drv{i})~=[nr{i},nv]                                         % MMP, 11/21/2026 (was)
    elseif ~isequal(size(Drv{i}),[nr{i},nv])                                % MMP, 11/21/2026
        error(['error in size of Drv',int2str(i)])
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%% Group 3 %%%%%%%%%%%%%%%%%%%%%%%

for i=1:nK
    % Cv{i}, Cvd{i}: zero where undefined or empty, else size checked
    % (nv x nr{i}) in both dimensions                                       % MMP, 11/21/2026
    if lcv<i||isempty(Cv{i})
        Cv{i}=zeros(nv,nr{i});
        disp(['Cv{',int2str(i),'} undefined, defaulting to 0'])
%   elseif size(Cv{i})~=[nv,nr{i}]                                          % MMP, 11/21/2026 (was)
    elseif ~isequal(size(Cv{i}),[nv,nr{i}])                                 % MMP, 11/21/2026
        error(['error in size of Cv',int2str(i)])
    end
    if lcvd<i||isempty(Cvd{i})
        Cvd{i}=zeros(nv,nr{i});
        disp(['Cvd{',int2str(i),'} undefined, defaulting to 0'])
%   elseif size(Cvd{i})~=[nv,nr{i}]                                         % MMP, 11/21/2026 (was)
    elseif ~isequal(size(Cvd{i}),[nv,nr{i}])                                % MMP, 11/21/2026
        error(['error in size of Cvd',int2str(i)])
    end
end


%%%%%%%%%%%%%%%%%%%%%%% Group 4 %%%%%%%%%%%%%%%%%%%%%%%


% C1: zero if undefined, else size checked (nz x nx) in both dimensions     % MMP, 11/21/2026
if ~exist('C1','var')
    C1=zeros(nz,nx);
    disp('C1 undefined, defaulting to 0')
% elseif size(C1)~=[nz,nx]                                                  % MMP, 11/21/2026 (was)
elseif ~isequal(size(C1),[nz,nx])                                           % MMP, 11/21/2026
    error(['error in size of C1'])
end

% D11: likewise (nz x nw)
if ~exist('D11','var')
    D11=zeros(nz,nw);
    disp('D11 undefined, defaulting to 0')
% elseif size(D11)~=[nz,nw]                                                 % MMP, 11/21/2026 (was)
elseif ~isequal(size(D11),[nz,nw])                                          % MMP, 11/21/2026
    error(['error in size of D11'])
end

% D12: likewise (nz x nu)
if ~exist('D12','var')
    D12=zeros(nz,nu);
    disp('D12 undefined, defaulting to 0')
% elseif size(D12)~=[nz,nu]                                                 % MMP, 11/21/2026 (was)
elseif ~isequal(size(D12),[nz,nu])                                          % MMP, 11/21/2026
    error(['error in size of D12'])
end

% D1v: likewise (nz x nv)
if ~exist('D1v','var')
    D1v=zeros(nz,nv);
    disp('D1v undefined, defaulting to 0')
% elseif size(D1v)~=[nz,nv]                                                 % MMP, 11/21/2026 (was)
elseif ~isequal(size(D1v),[nz,nv])                                          % MMP, 11/21/2026
    error(['error in size of D1v'])
end




%%%%%%%%%%%%%%%%%%%%%%% Group 5 %%%%%%%%%%%%%%%%%%%%%%%


% C2: zero if undefined, else size checked (ny x nx) in both dimensions     % MMP, 11/21/2026
if ~exist('C2','var')
    C2=zeros(ny,nx);
    disp('C2 undefined, defaulting to 0')
% elseif size(C2)~=[ny,nx]                                                  % MMP, 11/21/2026 (was)
elseif ~isequal(size(C2),[ny,nx])                                           % MMP, 11/21/2026
    error(['error in size of C2'])
end

% D21: likewise (ny x nw)
if ~exist('D21','var')
    D21=zeros(ny,nw);
    disp('D21 undefined, defaulting to 0')
% elseif size(D21)~=[ny,nw]                                                 % MMP, 11/21/2026 (was)
elseif ~isequal(size(D21),[ny,nw])                                          % MMP, 11/21/2026
    error(['error in size of D21'])
end

% D22: likewise (ny x nu)
if ~exist('D22','var')
    D22=zeros(ny,nu);
    disp('D22 undefined, defaulting to 0')
% elseif size(D22)~=[ny,nu]                                                 % MMP, 11/21/2026 (was)
elseif ~isequal(size(D22),[ny,nu])                                          % MMP, 11/21/2026
    error(['error in size of D22'])
end

% D2v: likewise (ny x nv)
if ~exist('D2v','var')
    D2v=zeros(ny,nv);
    disp('D2v undefined, defaulting to 0')
% elseif size(D2v)~=[ny,nv]                                                 % MMP, 11/21/2026 (was)
elseif ~isequal(size(D2v),[ny,nv])                                          % MMP, 11/21/2026
    error(['error in size of D2v'])
end



%clear lcr lbr1 lbr2 ldrv
%clear lcv lcvd

    % Output: the delays and the coefficients of x, z, y
    DDF_out.tau=tau;
    DDF_out.A0=A0;
    DDF_out.B1=B1;
    DDF_out.B2=B2;
    DDF_out.C1=C1;
    DDF_out.C2=C2;
    DDF_out.D11=D11;
    DDF_out.D12=D12;
    DDF_out.D21=D21;
    DDF_out.D22=D22;
    % the delay channels r_i and their outputs v
    DDF_out.Cr=Cr;
    DDF_out.Br1=Br1;
    DDF_out.Br2=Br2;
    DDF_out.Drv=Drv;
    DDF_out.Cv=Cv;
    DDF_out.Cvd=Cvd;
    DDF_out.Bv=Bv;
    DDF_out.D1v=D1v;
    DDF_out.D2v=D2v;




end


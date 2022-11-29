%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_DDF.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DDF_out=initialize_PIETOOLS_DDF(DDF)
% This a initialize/preprocessing routine, executed prior to 
% convert_PIETOOLS_DDF.m, which then converts a DDF to a PIE.
% The preprocessing stage checks which matrices have been defined and
% defines associated auxiliary variables used by the script
% convert_PIETOOLS_ddf.
% undefined signals are set to length 1 with 0 value
%
% input: DDF - a partially-defined DDF data structure
%
% output: DDF - a fully-defined DDF data structure
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
% Initial coding MP - 10_01_2020
% MMP, 5_30_2021; added new data structure and made script a function


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
%%%% arrays to make sure they are less than nK
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
        if exist('A0','var')
            nx=size(A0,1);
        elseif exist('B1','var')
            nx=size(B1,1,'var');
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

% In the case where no states are detected, make sure no states are used
% The places states are used is in C1, C2, and Cr{i}

% scanning for the number of regulated outputs
if ~exist('nz','var')
    if exist('C1','var')||exist('D11','var')||exist('D12','var')||exist('D1v','var')
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

if lcvd==0&&lcv==0
        warning(['Problem: The delayed channels have no outputs'])
        nv=0;
end    
nv=0;
for i=1:nK
    if lcv>=i&&~isempty(Cv{i})
        nv=size(Cv{i},1);
    elseif lcvd>=i&&~isempty(Cvd{i})
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
        disp('No disturbance channels detected')
        now=1;
        nw=0;
    end
end

% scanning for the number of control inputs
if ~exist('nu','var')
    if exist('B2','var')||exist('D12','var')||exist('D22','var')||exist('Br2','var')
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
if ~exist('A0','var')
    A0=zeros(nx,nx);
    disp('A0 undefined, defaulting to 0')
elseif (size(A0,1)~=nx)||(size(A0,2)~=nx)
    error('The dimensions of A0 are inconsistent with the detected number of inputs or states')
end


if ~exist('B1','var')
    B1=zeros(nx,nw);
    disp('B1 undefined, defaulting to 0')
elseif (size(B1,1)~=nx)||(size(B1,2)~=nw)
    error('The dimensions of B1 are inconsistent with the detected number of inputs or states')
end

if ~exist('B2','var')
    B2=zeros(nx,nu);
    disp('B2 undefined, defaulting to 0')
elseif (size(B2,1)~=nx)||(size(B2,2)~=nu)
    error('The dimensions of B2 are inconsistent with the detected number of inputs or states')
end

if ~exist('Bv','var')
    Bv=zeros(nx,nv);
    disp('Bv undefined, defaulting to 0')
elseif (size(Bv,1)~=nx)||(size(Bv,2)~=nv)
    error('The dimensions of Bv are inconsistent with the detected number of inputs or states')
end



%%%%%%%%%%%%%%%%%%%%%%% Group 2 %%%%%%%%%%%%%%%%%%%%%%%

for i=1:nK
    if lcr<i||isempty(Cr{i})
        Cr{i}=zeros(nr{i},nx);
        disp(['Cr{',int2str(i),'} undefined, defaulting to 0'])
    elseif size(Cr{i})~=[nr{i},nx]
        error(['error in size of Cr',int2str(i)])
    end
    
    if lbr1<i||isempty(Br1{i})
        Br1{i}=zeros(nr{i},nw);
        disp(['Br1{',int2str(i),'} undefined, defaulting to 0'])
    elseif size(Br1{i})~=[nr{i},nw]
        error(['error in size of Br1',int2str(i)])
    end
    
    if lbr2<i||isempty(Br2{i})
        Br2{i}=zeros(nr{i},nu);
        disp(['Br2{',int2str(i),'} undefined, defaulting to 0'])
    elseif size(Br2{i})~=[nr{i},nu]
        error(['error in size of Br2',int2str(i)])
    end
    
    if ldrv<i||isempty(Drv{i})
        Drv{i}=zeros(nr{i},nv);
        disp(['Drv{',int2str(i),'} undefined, defaulting to 0'])
    elseif size(Drv{i})~=[nr{i},nv]
        error(['error in size of Drv',int2str(i)])
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%% Group 3 %%%%%%%%%%%%%%%%%%%%%%%

for i=1:nK
    if lcv<i||isempty(Cv{i})
        Cv{i}=zeros(nv,nr{i});
        disp(['Cv{',int2str(i),'} undefined, defaulting to 0'])
    elseif size(Cv{i})~=[nv,nr{i}]
        error(['error in size of Cv',int2str(i)])
    end
    if lcvd<i||isempty(Cvd{i})
        Cvd{i}=zeros(nv,nr{i});
        disp(['Cvd{',int2str(i),'} undefined, defaulting to 0'])
    elseif size(Cvd{i})~=[nv,nr{i}]
        error(['error in size of Cvd',int2str(i)])
    end
end


%%%%%%%%%%%%%%%%%%%%%%% Group 4 %%%%%%%%%%%%%%%%%%%%%%%


if ~exist('C1','var')
    C1=zeros(nz,nx);
    disp('C1 undefined, defaulting to 0')
elseif size(C1)~=[nz,nx]
    error(['error in size of C1'])
end

if ~exist('D11','var')
    D11=zeros(nz,nw);
    disp('D11 undefined, defaulting to 0')
elseif size(D11)~=[nz,nw]
    error(['error in size of D11'])
end

if ~exist('D12','var')
    D12=zeros(nz,nu);
    disp('D12 undefined, defaulting to 0')
elseif size(D12)~=[nz,nu]
    error(['error in size of D12'])
end

if ~exist('D1v','var')
    D1v=zeros(nz,nv);
    disp('D1v undefined, defaulting to 0')
elseif size(D1v)~=[nz,nv]
    error(['error in size of D1v'])
end




%%%%%%%%%%%%%%%%%%%%%%% Group 5 %%%%%%%%%%%%%%%%%%%%%%%


if ~exist('C2','var')
    C2=zeros(ny,nx);
    disp('C2 undefined, defaulting to 0')
elseif size(C2)~=[ny,nx]
    error(['error in size of C2'])
end

if ~exist('D21','var')
    D21=zeros(ny,nw);
    disp('D21 undefined, defaulting to 0')
elseif size(D21)~=[ny,nw]
    error(['error in size of D21'])
end

if ~exist('D22','var')
    D22=zeros(ny,nu);
    disp('D22 undefined, defaulting to 0')
elseif size(D22)~=[ny,nu]
    error(['error in size of D22'])
end

if ~exist('D2v','var')
    D2v=zeros(ny,nv);
    disp('D2v undefined, defaulting to 0')
elseif size(D2v)~=[ny,nv]
    error(['error in size of D2v'])
end



%clear lcr lbr1 lbr2 ldrv
%clear lcv lcvd

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

% Now for the hard part - setting up the PIE

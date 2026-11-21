%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_NDS.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NDS_out=initialize_PIETOOLS_NDS(NDS)
% initializes and checks the NDS formulation, with K = numel(tau) delays,
%   \dot x(t) = A0 x(t) + sum_i Ai{i} x(t-tau(i)) + sum_i int_{-tau(i)}^0 Adi{i}(s) x(t+s) ds
%             + sum_i Ei{i} \dot x(t-tau(i)) + sum_i int_{-tau(i)}^0 Edi{i}(s) \dot x(t+s) ds
%             + B1 w(t) + sum_i B1i{i} w(t-tau(i)) + sum_i int_{-tau(i)}^0 B1di{i}(s) w(t+s) ds
%             + B2 u(t) + sum_i B2i{i} u(t-tau(i)) + sum_i int_{-tau(i)}^0 B2di{i}(s) u(t+s) ds
%        z(t) = C1 x(t) + D11 w(t) + D12 u(t) + (C1i, C1di, C1ei, C1dei, D11i, ..., D12di likewise)
%        y(t) = C2 x(t) + D21 w(t) + D22 u(t) + (C2i, C2di, C2ei, C2dei, D21i, ..., D22di likewise)
% where C1ei, C1dei (C2ei, C2dei) multiply \dot x(t-tau(i)) and \dot x(t+s).
% undefined signals are set to length 0 with 0 value.
% The dimensions nx, nw, nu, nz, ny are read off the first defined
% coefficient of each signal (of a cell, its last element); undefined
% coefficients are set to zero, defined ones are checked against them.
%
% input: NDS - a partially-defined NDS data structure: 'struct' with the
%        field tau (1 x K 'double', delays >= 0; required), the 'double'
%        matrices A0, B1, B2, C1, C2, D11, D12, D21, D22, cells of at
%        most K 'double' matrices Ai, Ei, B1i, ..., C1ei, C2ei, ..., D22i
%        and of 'double' or 'polynomial'
%        kernels Adi, Edi, B1di, ..., C1dei, C2dei, ..., D22di (empty or
%        missing elements are zero), and optional fields nx, nw, nu, nz,
%        ny that fix the dimensions.
%
% output: NDS_out - a fully-defined NDS data structure: 'struct' with the
%        fields tau, A0, ..., D22 and the 1 x K cells Ai, ..., D22di, Ei,
%        Edi, C1ei, C2ei, C1dei, C2dei, each of its size.

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
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/20/2026: A0 defaulted and checked; dimensions read off Ei, Edi, C1ei,
%                   C1dei, C2ei, C2dei, B1di (unreached); field names in errors.


%%%%%%%%%%%%%%%%%%%%%%%%%
% Dumping the data structure to the workspace. We are not bad people!
names=fieldnames(NDS);
for i=1:numel(names)
    eval([names{i},'=NDS.',names{i},';'])
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% checking to make sure the delays are defined.
if ~exist('tau','var')
    error('Ooops, you forgot to define the delays in tau. Otherwise, what is the point? Go use an ODE solver.')
end
nK=length(tau);    % the specified number of delays

for i=1:nK
    if tau(i)<0
        error(['tau(',int2str(i),') is negative. Delay values must be positive.'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Before Proceeding, we should check that the user has defined enough 
%%%% delays. For this we check the number of elements in each of the cell
%%%% arrays to make sure they are at most nK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State equation: Ai, Adi, B1i, B1di, B2i, B2di
if exist('Ai','var')&&(numel(Ai)>nK)
    error('You have defined an element of Ai for which there is no matching value of delay.')
end

if exist('Adi','var')&&(numel(Adi)>nK)
    error('You have defined an element of Adi for which there is no matching value of delay.')
end

if exist('B1i','var')&&(numel(B1i)>nK)
    error('You have defined an element of B1i for which there is no matching value of delay.')
end

if exist('B1di','var')&&(numel(B1di)>nK)
    error('You have defined an element of B1di for which there is no matching value of delay.')
end

if exist('B2i','var')&&(numel(B2i)>nK)
    error('You have defined an element of B2i for which there is no matching value of delay.')
end

if exist('B2di','var')&&(numel(B2di)>nK)
    error('You have defined an element of B2di for which there is no matching value of delay.')
end

% Output equations, state terms: C1i, C1di, C2i, C2di
if exist('C1i','var')&&(numel(C1i)>nK)
    error('You have defined an element of C1i for which there is no matching value of delay.')
end
if exist('C1di','var')&&(numel(C1di)>nK)
    error('You have defined an element of C1di for which there is no matching value of delay.')
end

if exist('C2i','var')&&(numel(C2i)>nK)
    error('You have defined an element of C2i for which there is no matching value of delay.')
end
if exist('C2di','var')&&(numel(C2di)>nK)
    error('You have defined an element of C2di for which there is no matching value of delay.')
end

% Regulated output, input terms: D11i, D11di, D12i, D12di
if exist('D11i','var')&&(numel(D11i)>nK)
    error('You have defined an element of D11i for which there is no matching value of delay.')
end

if exist('D11di','var')&&(numel(D11di)>nK)
    error('You have defined an element of D11di for which there is no matching value of delay.')
end

if exist('D12i','var')&&(numel(D12i)>nK)
    error('You have defined an element of D12i for which there is no matching value of delay.')
end

if exist('D12di','var')&&(numel(D12di)>nK)
    error('You have defined an element of D12di for which there is no matching value of delay.')
end

% Sensed output, input terms: D21i, D21di, D22i, D22di
if exist('D21i','var')&&(numel(D21i)>nK)
    error('You have defined an element of D21i for which there is no matching value of delay.')
end

if exist('D21di','var')&&(numel(D21di)>nK)
    error('You have defined an element of D21di for which there is no matching value of delay.')
end

if exist('D22i','var')&&(numel(D22i)>nK)
    error('You have defined an element of D22i for which there is no matching value of delay.')
end

if exist('D22di','var')&&(numel(D22di)>nK)
    error('You have defined an element of D22di for which there is no matching value of delay.')
end

%%%%%%%%%%% New for NDS %%%%%%%%%%%%%%%%%%%%%%
if exist('Ei','var')&&(numel(Ei)>nK)
    error('You have defined an element of Ei for which there is no matching value of delay.')
end

if exist('C1ei','var')&&(numel(C1ei)>nK)
%   error('You have defined an element of Ei for which there is no matching value of delay.') % MMP, 11/20/2026 (was)
    error('You have defined an element of C1ei for which there is no matching value of delay.') % MMP, 11/20/2026
end

if exist('C2ei','var')&&(numel(C2ei)>nK)
%   error('You have defined an element of Ei for which there is no matching value of delay.') % MMP, 11/20/2026 (was)
    error('You have defined an element of C2ei for which there is no matching value of delay.') % MMP, 11/20/2026
end


if exist('Edi','var')&&(numel(Edi)>nK)
%   error('You have defined an element of Ei for which there is no matching value of delay.') % MMP, 11/20/2026 (was)
    error('You have defined an element of Edi for which there is no matching value of delay.') % MMP, 11/20/2026
end

if exist('C1dei','var')&&(numel(C1dei)>nK)
%   error('You have defined an element of Ei for which there is no matching value of delay.') % MMP, 11/20/2026 (was)
    error('You have defined an element of C1dei for which there is no matching value of delay.') % MMP, 11/20/2026
end

if exist('C2dei','var')&&(numel(C2dei)>nK)
%   error('You have defined an element of Ei for which there is no matching value of delay.') % MMP, 11/20/2026 (was)
    error('You have defined an element of C2dei for which there is no matching value of delay.') % MMP, 11/20/2026
end
%%%%%%%%%%% End New for NDS %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Here we detect the number of states, assuming the user has
%%%% not defined these inconsistently. If they have defined them
%%%% inconsistently, this will generate a hard-to-trace error in the
%%%% conversion script. The error checking is fairly robust, however.
%%%% If only the cellular input are defined, it does assume the last
%%%% cellular entry is defined by the user.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nstates=0;noregs=0;nomeas=0;now=0;nou=0; % toggles to determine if no states, inputs or outputs are present

% scanning for the number of states
if ~exist('nx','var')
    % Ei, Edi included, else their branches below are never reached         % MMP, 11/20/2026
%   if exist('A0','var')||exist('Ai','var')||exist('Adi','var')||exist('B1','var')||exist('B1i','var')||exist('B1di','var')||exist('B2','var')||exist('B2i','var')||exist('B2di','var') % MMP, 11/20/2026 (was)
    if exist('A0','var')||exist('Ai','var')||exist('Adi','var')||exist('B1','var')||exist('B1i','var')||exist('B1di','var')||exist('B2','var')||exist('B2i','var')||exist('B2di','var')||exist('Ei','var')||exist('Edi','var') % MMP, 11/20/2026
        % nx: rows of A0, B1 or B2, the first one defined
        if exist('A0','var')
            nx=size(A0,1);
        elseif exist('B1','var')
            nx=size(B1,1);
        elseif exist('B2','var')
            nx=size(B2,1);
        % else rows of the last element of Ai, B1i, B2i, Adi, B1di, B2di,   % MMP, 11/20/2026
        % Ei or Edi                                                         % MMP, 11/20/2026
        elseif exist('Ai','var')
            nx=size(Ai{numel(Ai)},1);
        elseif exist('B1i','var')
            nx=size(B1i{numel(B1i)},1);
        elseif exist('B2i','var')
            nx=size(B2i{numel(B2i)},1);
        elseif exist('Adi','var')
            nx=size(Adi{numel(Adi)},1);
        elseif exist('B1di','var')
            nx=size(B1di{numel(B1di)},1);
        elseif exist('B2di','var')
            nx=size(B2di{numel(B2di)},1);
        elseif exist('Ei','var')
            nx=size(Ei{numel(Ei)},1);
        elseif exist('Edi','var')
            nx=size(Edi{numel(Edi)},1);
        end
    else
        % No state-equation coefficient: no states
        nx=0;
        warning('No state dynamics detected. There are no states')
        nstates=1;
    end
end

% A0: zero if undefined, else size checked (nx x nx), as in                 % MMP, 11/20/2026
% initialize_PIETOOLS_DDE; without it NDS_out.A0=A0 failed on a missing A0  % MMP, 11/20/2026
if ~exist('A0','var')                                                       % MMP, 11/20/2026
    A0=zeros(nx,nx);                                                        % MMP, 11/20/2026
    disp('A0 undefined, defaulting to 0')                                   % MMP, 11/20/2026
elseif (size(A0,1)~=nx)||(size(A0,2)~=nx)                                   % MMP, 11/20/2026
    error('The dimensions of A0 are inconsistent with the detected number of states') % MMP, 11/20/2026
end                                                                         % MMP, 11/20/2026


% Ai{i}, i = 1..nK: zero where undefined or empty, else size checked (nx x nx)
if ~exist('Ai','var')
    for i=1:nK
        Ai{i}=zeros(nx);
    end
    disp('Ai undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(Ai)<i)||isempty(Ai{i})
            Ai{i}=zeros(nx);
        elseif (size(Ai{i},1)~=nx)||(size(Ai{i},2)~=nx)
            error(['The dimensions of Ai{',int2str(i),'} are inconsistent with the detected number of states'])
        end
    end
end

% Adi{i}: zero where undefined or empty, else size checked (nx x nx)
if ~exist('Adi','var')
    for i=1:nK
        Adi{i}=zeros(nx);
    end
    disp('Adi undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(Adi)<i)||isempty(Adi{i})
            Adi{i}=zeros(nx);
        elseif (size(Adi{i},1)~=nx)||(size(Adi{i},2)~=nx)
            error(['The dimensions of Adi{',int2str(i),'} are inconsistent with the detected number of states'])
        end
    end
end


%%%%%%%%%%% New for NDS %%%%%%%%%%%%%%%%%%%%%%
% Ei{i}: zero where undefined or empty, else size checked (nx x nx)
if ~exist('Ei','var')
    for i=1:nK
        Ei{i}=zeros(nx);
    end
    disp('Ei undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(Ei)<i)||isempty(Ei{i})
            Ei{i}=zeros(nx);
        elseif (size(Ei{i},1)~=nx)||(size(Ei{i},2)~=nx)
            error(['The dimensions of Ei{',int2str(i),'} are inconsistent with the detected number of states'])
        end
    end
end

% Edi{i}: zero where undefined or empty, else size checked (nx x nx)
if ~exist('Edi','var')
    for i=1:nK
        Edi{i}=zeros(nx);
    end
    disp('Edi undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(Edi)<i)||isempty(Edi{i})
            Edi{i}=zeros(nx);
        elseif (size(Edi{i},1)~=nx)||(size(Edi{i},2)~=nx)
            error(['The dimensions of Edi{',int2str(i),'} are inconsistent with the detected number of states'])
        end
    end
end



% if ~exist('Ai')
%     for i=1:nK
%         Ai{i}=zeros(nx);
%     end
%     disp('No discrete delays detected. Defaulting to 0')
% end
% 
% 
% if numel(Ai)<nK
%     for i=(numel(Ai)+1):nK
%         Ai{i}=zeros(nx);
%         disp('Missing discrete delay term. Defaulting to 0')
%     end
% end
% 
% for i=1:nK
%     if size(Ai{i})~=size(A0)
%         Ai{i}=zeros(nx);
%         disp('Missing or malformed discrete delay term. Resetting to 0')
%     end
% end
% 
% if ~exist('Ad')
%     for i=1:nK
%         Ad{i}=zeros(nx);
%     end
%     disp('No distributed delays detecting. Defaulting to 0')
% end
%     
% 
% for i=1:nK
%     if size(Ad{i})~=size(A0)
%         Ad{i}=zeros(nx);
%         disp('Missing or malformed distributed delay term. Resetting to 0')
%     end
% end



% if ~exist('B1')
%     B1=zeros(nx,1);
%     disp('B1 undefined, defaulting to 0')
% end
% nw=size(B1,2);   % number of disturbances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Here we detect the number of inputs and outputs, assuming the user has
%%%% not defined these inconsistently. If they have defined them
%%%% inconsistently, this will generate a hard-to-trace error in the
%%%% conversion script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% scanning for the number of disturbance inputs
if ~exist('nw','var')
    % B1di, the field name (Bdi is none)                                    % MMP, 11/20/2026
%   if exist('B1','var')||exist('B1i','var')||exist('Bdi','var')||exist('D11','var')||exist('D11i','var')||exist('D11di','var')||exist('D21','var')||exist('D21i','var')||exist('D21di','var') % MMP, 11/20/2026 (was)
    if exist('B1','var')||exist('B1i','var')||exist('B1di','var')||exist('D11','var')||exist('D11i','var')||exist('D11di','var')||exist('D21','var')||exist('D21i','var')||exist('D21di','var') % MMP, 11/20/2026
        % nw: columns of B1, D11 or D21, the first one defined
        if exist('B1')
            nw=size(B1,2);
        elseif exist('D11','var')
            nw=size(D11,2);
        elseif exist('D21','var')
            nw=size(D21,2);
        % else columns of the last element of B1i, D11i, D21i, B1di, D11di or D21di
        elseif exist('B1i','var')
            nw=size(B1i{numel(B1i)},2);
        elseif exist('D11i','var')
            nw=size(D11i{numel(D11i)},2);
        elseif exist('D21i','var')
            nw=size(D21i{numel(D21i)},2);
        elseif exist('B1di','var')
            nw=size(B1di{numel(B1di)},2);
        elseif exist('D11di','var')
            nw=size(D11di{numel(D11di)},2);
        elseif exist('D21di','var')
            nw=size(D21di{numel(D21di)},2);
        end
    else
        % No disturbance coefficient: nw = 0
        now=1;
        nw=0;
    end
end


% scanning for the number of control inputs
if ~exist('nu','var')
    if exist('B2','var')||exist('B2i','var')||exist('B2di','var')||exist('D12','var')||exist('D12i','var')||exist('D12di','var')||exist('D22','var')||exist('D22i','var')||exist('D22di','var')
        % nu: columns of B2, D12 or D22, the first one defined
        if exist('B2','var')
            nu=size(B2,2);
        elseif exist('D12','var')
            nu=size(D12,2);
        elseif exist('D22','var')
            nu=size(D22,2);
        % else columns of the last element of B2i, D12i, D22i, B2di, D12di or D22di
        elseif exist('B2i','var')
            nu=size(B2i{numel(B2i)},2);
        elseif exist('D12i','var')
            nu=size(D12i{numel(D12i)},2);
        elseif exist('D22i','var')
            nu=size(D22i{numel(D22i)},2);
        elseif exist('B2di','var')
            nu=size(B2di{numel(B2di)},2);
        elseif exist('D12di','var')
            nu=size(D12di{numel(D12di)},2);
        elseif exist('D22di','var')
            nu=size(D22di{numel(D22di)},2);
        end
    else
        % No control-input coefficient: nu = 0
        nu=0;
        nou=1;
    end
end

% scanning for the number of regulated outputs
if ~exist('nz','var')
    % C1ei, C1dei included, else their branches below are never reached     % MMP, 11/20/2026
%   if exist('C1','var')||exist('C1i','var')||exist('C1di','var')||exist('D11','var')||exist('D11i','var')||exist('D11di','var')||exist('D12','var')||exist('D12i','var')||exist('D12di','var') % MMP, 11/20/2026 (was)
    if exist('C1','var')||exist('C1i','var')||exist('C1di','var')||exist('D11','var')||exist('D11i','var')||exist('D11di','var')||exist('D12','var')||exist('D12i','var')||exist('D12di','var')||exist('C1ei','var')||exist('C1dei','var') % MMP, 11/20/2026
        % nz: rows of C1, D11 or D12, the first one defined
        if exist('C1','var')
            nz=size(C1,1);
        elseif exist('D11','var')
            nz=size(D11,1);
        elseif exist('D12','var')
            nz=size(D12,1);
        % else rows of the last element of C1i, D11i, D12i, C1di, D11di,    % MMP, 11/20/2026
        % D12di, C1ei or C1dei                                              % MMP, 11/20/2026
        elseif exist('C1i','var')
            nz=size(C1i{numel(C1i)},1);
        elseif exist('D11i','var')
            nz=size(D11i{numel(D11i)},1);
        elseif exist('D12i','var')
            nz=size(D12i{numel(D12i)},1);
        elseif exist('C1di','var')
            nz=size(C1di{numel(C1di)},1);
        elseif exist('D11di','var')
            nz=size(D11di{numel(D11di)},1);
        elseif exist('D12di','var')
            nz=size(D12di{numel(D12di)},1);
        elseif exist('C1ei','var')
            nz=size(C1ei{numel(C1ei)},1);
        elseif exist('C1dei','var')
            nz=size(C1dei{numel(C1dei)},1);
        end
    else
        % No regulated-output coefficient: nz = 0
        nz=0;
        noregs=1;
    end
end


% scanning for the number of sensed outputs
if ~exist('ny','var')
    % C2ei, C2dei included, else their branches below are never reached     % MMP, 11/20/2026
%   if exist('C2','var')||exist('C2i','var')||exist('C2di','var')||exist('D21','var')||exist('D21i','var')||exist('D21di','var')||exist('D22','var')||exist('D22i','var')||exist('D22di','var') % MMP, 11/20/2026 (was)
    if exist('C2','var')||exist('C2i','var')||exist('C2di','var')||exist('D21','var')||exist('D21i','var')||exist('D21di','var')||exist('D22','var')||exist('D22i','var')||exist('D22di','var')||exist('C2ei','var')||exist('C2dei','var') % MMP, 11/20/2026
        % ny: rows of C2, D21 or D22, the first one defined
        if exist('C2','var')
            ny=size(C2,1);
        elseif exist('D21','var')
            ny=size(D21,1);
        elseif exist('D22','var')
            ny=size(D22,1);
        % else rows of the last element of C2i, D21i, D22i, C2di, D21di,    % MMP, 11/20/2026
        % D22di, C2ei or C2dei                                              % MMP, 11/20/2026
        elseif exist('C2i','var')
            ny=size(C2i{numel(C2i)},1);
        elseif exist('D21i','var')
            ny=size(D21i{numel(D21i)},1);
        elseif exist('D22i','var')
            ny=size(D22i{numel(D22i)},1);
        elseif exist('C2di','var')
            ny=size(C2di{numel(C2di)},1);
        elseif exist('D21di','var')
            ny=size(D21di{numel(D21di)},1);
        elseif exist('D22di','var')
            ny=size(D22di{numel(D22di)},1);
        elseif exist('C2ei','var')
%           nz=size(C2ei{numel(C2ei)},1);                                   % MMP, 11/20/2026 (was)
            ny=size(C2ei{numel(C2ei)},1);                                   % MMP, 11/20/2026
        elseif exist('C2dei','var')
%           nz=size(C2dei{numel(C2dei)},1);                                 % MMP, 11/20/2026 (was)
            ny=size(C2dei{numel(C2dei)},1);                                 % MMP, 11/20/2026
        end
    else
        % No sensed-output coefficient: ny = 0
        ny=0;
        nomeas=1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Now that we know the number of inputs and outputs, we can check each
%%%% of the user-defined matrices. If they don't exist they are defaulted
%%%% to zero. If they exist, we verify they are of the correct dimension.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I/O terms in the dynamics
% B1: zero if undefined, else size checked (nx x nw)
if ~exist('B1','var')
    B1=zeros(nx,nw);
    disp('B1 undefined, defaulting to 0')
elseif (size(B1,1)~=nx)||(size(B1,2)~=nw)
    error('The dimensions of B1 are inconsistent with the detected number of inputs or states')
end

% B1i{i}: zero where undefined or empty, else size checked (nx x nw)
if ~exist('B1i','var')
    for i=1:nK
        B1i{i}=zeros(nx,nw);
    end
    disp('B1i undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(B1i)<i)||isempty(B1i{i})
            B1i{i}=zeros(nx,nw);
        elseif (size(B1i{i},1)~=nx)||(size(B1i{i},2)~=nw)
            error(['The dimensions of B1i{',int2str(i),'} are inconsistent with the detected number of inputs or states'])
        end
    end
end

% B1di{i}: zero where undefined or empty, else size checked (nx x nw)
if ~exist('B1di','var')
    for i=1:nK
        B1di{i}=zeros(nx,nw);
    end
    disp('B1di undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(B1di)<i)||isempty(B1di{i})
            B1di{i}=zeros(nx,nw);
        elseif (size(B1di{i},1)~=nx)||(size(B1di{i},2)~=nw)
            error(['The dimensions of B1di{',int2str(i),'} are inconsistent with the detected number of inputs or states'])
        end
    end
end

% B2: zero if undefined, else size checked (nx x nu)
if ~exist('B2','var')
    B2=zeros(nx,nu);
    disp('B2 undefined, defaulting to 0')
elseif (size(B2,1)~=nx)||(size(B2,2)~=nu)
    error('The dimensions of B2 are inconsistent with the detected number of inputs or states')
end

% B2i{i}: zero where undefined or empty, else size checked (nx x nu)
if ~exist('B2i','var')
    for i=1:nK
        B2i{i}=zeros(nx,nu);
    end
    disp('B2i undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(B2i)<i)||isempty(B2i{i})
            B2i{i}=zeros(nx,nu);
        elseif (size(B2i{i},1)~=nx)||(size(B2i{i},2)~=nu)
            error(['The dimensions of B2i{',int2str(i),'} are inconsistent with the detected number of inputs or states'])
        end
    end
end

% B2di{i}: zero where undefined or empty, else size checked (nx x nu)
if ~exist('B2di','var')
    for i=1:nK
        B2di{i}=zeros(nx,nu);
    end
    disp('B2di undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(B2di)<i)||isempty(B2di{i})
            B2di{i}=zeros(nx,nu);
        elseif (size(B2di{i},1)~=nx)||(size(B2di{i},2)~=nu)
            error(['The dimensions of B2di{',int2str(i),'} are inconsistent with the detected number of inputs or states'])
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Now error checking for the regulated output expression
% C1: zero if undefined, else size checked (nz x nx)
if ~exist('C1','var')
    C1=zeros(nz,nx);
    disp('C1 undefined, defaulting to 0')
elseif (size(C1,1)~=nz)||(size(C1,2)~=nx)
    error('The dimensions of C1 are inconsistent with the detected number of regulated outputs or states')
end

% C1i{i}: zero where undefined or empty, else size checked (nz x nx)
if ~exist('C1i','var')
    for i=1:nK
        C1i{i}=zeros(nz,nx);
    end
    disp('C1i undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(C1i)<i)||isempty(C1i{i})
            C1i{i}=zeros(nz,nx);
        elseif (size(C1i{i},1)~=nz)||(size(C1i{i},2)~=nx)
            error(['The dimensions of C1i{',int2str(i),'} are inconsistent with the detected number of regulated outputs or states'])
        end
    end
end

% C1di{i}: zero where undefined or empty, else size checked (nz x nx)
if ~exist('C1di','var')
    for i=1:nK
        C1di{i}=zeros(nz,nx);
    end
    disp('C1di undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(C1di)<i)||isempty(C1di{i})
            C1di{i}=zeros(nz,nx);
        elseif (size(C1di{i},1)~=nz)||(size(C1di{i},2)~=nx)
            error(['The dimensions of C1di{',int2str(i),'} are inconsistent with the detected number of regulated outputs or states'])
        end
    end
end

% C1ei{i}: zero where undefined or empty, else size checked (nz x nx)
if ~exist('C1ei','var')
    for i=1:nK
        C1ei{i}=zeros(nz,nx);
    end
    disp('C1ei undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(C1ei)<i)||isempty(C1ei{i})
            C1ei{i}=zeros(nz,nx);
        elseif (size(C1ei{i},1)~=nz)||(size(C1ei{i},2)~=nx)
            error(['The dimensions of C1ei{',int2str(i),'} are inconsistent with the detected number of regulated outputs or states'])
        end
    end
end

% C1dei{i}: zero where undefined or empty, else size checked (nz x nx)
if ~exist('C1dei','var')
    for i=1:nK
        C1dei{i}=zeros(nz,nx);
    end
    disp('C1dei undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(C1dei)<i)||isempty(C1dei{i})
            C1dei{i}=zeros(nz,nx);
        elseif (size(C1dei{i},1)~=nz)||(size(C1dei{i},2)~=nx)
            error(['The dimensions of C1dei{',int2str(i),'} are inconsistent with the detected number of regulated outputs or states'])
        end
    end
end





% D11: zero if undefined, else size checked (nz x nw)
if ~exist('D11','var')
    D11=zeros(nz,nw);
    disp('D11 undefined, defaulting to 0')
elseif (size(D11,1)~=nz)||(size(D11,2)~=nw)
    error('The dimensions of D11 are inconsistent with the detected number of regulated outputs or disturbances')
end

% D11i{i}: zero where undefined or empty, else size checked (nz x nw)
if ~exist('D11i','var')
    for i=1:nK
        D11i{i}=zeros(nz,nw);
    end
    disp('D11i undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(D11i)<i)||isempty(D11i{i})
            D11i{i}=zeros(nz,nw);
        elseif (size(D11i{i},1)~=nz)||(size(D11i{i},2)~=nw)
            error(['The dimensions of D11i{',int2str(i),'} are inconsistent with the detected number of regulated outputs or disturbances'])
        end
    end
end

% D11di{i}: zero where undefined or empty, else size checked (nz x nw)
if ~exist('D11di','var')
    for i=1:nK
        D11di{i}=zeros(nz,nw);
    end
    disp('D11di undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(D11di)<i)||isempty(D11di{i})
            D11di{i}=zeros(nz,nw);
        elseif (size(D11di{i},1)~=nz)||(size(D11di{i},2)~=nw)
            error(['The dimensions of D11di{',int2str(i),'} are inconsistent with the detected number of regulated outputs or disturbances'])
        end
    end
end

% D12: zero if undefined, else size checked (nz x nu)
if ~exist('D12','var')
    D12=zeros(nz,nu);
    disp('D12 undefined, defaulting to 0')
elseif (size(D12,1)~=nz)||(size(D12,2)~=nu)
    error('The dimensions of D12 are inconsistent with the detected number of regulated outputs or control inputs')
end

% D12i{i}: zero where undefined or empty, else size checked (nz x nu)
if ~exist('D12i','var')
    for i=1:nK
        D12i{i}=zeros(nz,nu);
    end
    disp('D12i undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(D12i)<i)||isempty(D12i{i})
            D12i{i}=zeros(nz,nu);
        elseif (size(D12i{i},1)~=nz)||(size(D12i{i},2)~=nu)
            error(['The dimensions of D12i{',int2str(i),'} are inconsistent with the detected number of regulated outputs or control inputs'])
        end
    end
end

% D12di{i}: zero where undefined or empty, else size checked (nz x nu)
if ~exist('D12di','var')
    for i=1:nK
        D12di{i}=zeros(nz,nu);
    end
    disp('D12di undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(D12di)<i)||isempty(D12di{i})
            D12di{i}=zeros(nz,nu);
        elseif (size(D12di{i},1)~=nz)||(size(D12di{i},2)~=nu)
            error(['The dimensions of D12di{',int2str(i),'} are inconsistent with the detected number of regulated outputs or control inputs'])
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Now error checking for the sensed output expression
% C2: zero if undefined, else size checked (ny x nx)
if ~exist('C2','var')
    C2=zeros(ny,nx);
    disp('C2 undefined, defaulting to 0')
elseif (size(C2,1)~=ny)||(size(C2,2)~=nx)
    error('The dimensions of C2 are inconsistent with the detected number of sensed outputs or states')
end

% C2i{i}: zero where undefined or empty, else size checked (ny x nx)
if ~exist('C2i','var')
    for i=1:nK
        C2i{i}=zeros(ny,nx);
    end
    disp('C2i undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(C2i)<i)||isempty(C2i{i})
            C2i{i}=zeros(ny,nx);
        elseif (size(C2i{i},1)~=ny)||(size(C2i{i},2)~=nx)
            error(['The dimensions of C2i{',int2str(i),'} are inconsistent with the detected number of sensed outputs or states'])
        end
    end
end

% C2di{i}: zero where undefined or empty, else size checked (ny x nx)
if ~exist('C2di','var')
    for i=1:nK
        C2di{i}=zeros(ny,nx);
    end
    disp('C2di undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(C2di)<i)||isempty(C2di{i})
            C2di{i}=zeros(ny,nx);
        elseif (size(C2di{i},1)~=ny)||(size(C2di{i},2)~=nx)
            error(['The dimensions of C2di{',int2str(i),'} are inconsistent with the detected number of sensed outputs or states'])
        end
    end
end

% C2ei{i}: zero where undefined or empty, else size checked (ny x nx)
if ~exist('C2ei','var')
    for i=1:nK
        C2ei{i}=zeros(ny,nx);
    end
    disp('C2ei undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(C2ei)<i)||isempty(C2ei{i})
            C2ei{i}=zeros(ny,nx);
        elseif (size(C2ei{i},1)~=ny)||(size(C2ei{i},2)~=nx)
            error(['The dimensions of C2ei{',int2str(i),'} are inconsistent with the detected number of sensed outputs or states'])
        end
    end
end

% C2dei{i}: zero where undefined or empty, else size checked (ny x nx)
if ~exist('C2dei','var')
    for i=1:nK
        C2dei{i}=zeros(ny,nx);
    end
    disp('C2dei undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(C2dei)<i)||isempty(C2dei{i})
            C2dei{i}=zeros(ny,nx);
        elseif (size(C2dei{i},1)~=ny)||(size(C2dei{i},2)~=nx)
            error(['The dimensions of C2dei{',int2str(i),'} are inconsistent with the detected number of sensed outputs or states'])
        end
    end
end




% D21: zero if undefined, else size checked (ny x nw)
if ~exist('D21','var')
    D21=zeros(ny,nw);
    disp('D21 undefined, defaulting to 0')
elseif (size(D21,1)~=ny)||(size(D21,2)~=nw)
    error('The dimensions of D21 are inconsistent with the detected number of sensed outputs or disturbances')
end

% D21i{i}: zero where undefined or empty, else size checked (ny x nw)
if ~exist('D21i','var')
    for i=1:nK
        D21i{i}=zeros(ny,nw);
    end
    disp('D21i undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(D21i)<i)||isempty(D21i{i})
            D21i{i}=zeros(ny,nw);
        elseif (size(D21i{i},1)~=ny)||(size(D21i{i},2)~=nw)
            error(['The dimensions of D21i{',int2str(i),'} are inconsistent with the detected number of sensed outputs or disturbances'])
        end
    end
end

% D21di{i}: zero where undefined or empty, else size checked (ny x nw)
if ~exist('D21di','var')
    for i=1:nK
        D21di{i}=zeros(ny,nw);
    end
    disp('D21di undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(D21di)<i)||isempty(D21di{i})
            D21di{i}=zeros(ny,nw);
        elseif (size(D21di{i},1)~=ny)||(size(D21di{i},2)~=nw)
            error(['The dimensions of D21di{',int2str(i),'} are inconsistent with the detected number of sensed outputs or disturbances'])
        end
    end
end

% D22: zero if undefined, else size checked (ny x nu)
if ~exist('D22','var')
    D22=zeros(ny,nu);
    disp('22 undefined, defaulting to 0')
elseif (size(D22,1)~=ny)||(size(D22,2)~=nu)
    error('The dimensions of D22 are inconsistent with the detected number of sensed outputs or control inputs')
end

% D22i{i}: zero where undefined or empty, else size checked (ny x nu)
if ~exist('D22i','var')
    for i=1:nK
        D22i{i}=zeros(ny,nu);
    end
    disp('D22i undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(D22i)<i)||isempty(D22i{i})
            D22i{i}=zeros(ny,nu);
        elseif (size(D22i{i},1)~=ny)||(size(D22i{i},2)~=nu)
            error(['The dimensions of D22i{',int2str(i),'} are inconsistent with the detected number of sensed outputs or control inputs'])
        end
    end
end

% D22di{i}: zero where undefined or empty, else size checked (ny x nu)
if ~exist('D22di','var')
    for i=1:nK
        D22di{i}=zeros(ny,nu);
    end
    disp('D22di undefined, defaulting to 0')
else
    for i=1:nK
        if (numel(D22di)<i)||isempty(D22di{i})
            D22di{i}=zeros(ny,nu);
        elseif (size(D22di{i},1)~=ny)||(size(D22di{i},2)~=nu)
            error(['The dimensions of D22di{',int2str(i),'} are inconsistent with the detected number of sensed outputs or control inputs'])
        end
    end
end

% Output: the delays and the coefficients of the terms at time t
NDS_out.tau=tau;
NDS_out.A0=A0;
NDS_out.B1=B1;
NDS_out.B2=B2;
NDS_out.C1=C1;
NDS_out.C2=C2;
NDS_out.D11=D11;
NDS_out.D12=D12;
NDS_out.D21=D21;
NDS_out.D22=D22;
    % the discrete delay coefficients
    NDS_out.Ai=Ai;
    NDS_out.B1i=B1i;
    NDS_out.B2i=B2i;
    NDS_out.C1i=C1i;
    NDS_out.C2i=C2i;
    NDS_out.D11i=D11i;
    NDS_out.D12i=D12i;
    NDS_out.D21i=D21i;
    NDS_out.D22i=D22i;
    % the distributed delay kernels
    NDS_out.Adi=Adi;
    NDS_out.B1di=B1di;
    NDS_out.B2di=B2di;
    NDS_out.C1di=C1di;
    NDS_out.C2di=C2di;
    NDS_out.D11di=D11di;
    NDS_out.D12di=D12di;
    NDS_out.D21di=D21di;
    NDS_out.D22di=D22di;
% new for NDS    
    NDS_out.Ei=Ei;
    NDS_out.Edi=Edi;
    NDS_out.C1ei=C1ei;
    NDS_out.C2ei=C2ei;
    NDS_out.C1dei=C1dei;
    NDS_out.C2dei=C2dei;

    end






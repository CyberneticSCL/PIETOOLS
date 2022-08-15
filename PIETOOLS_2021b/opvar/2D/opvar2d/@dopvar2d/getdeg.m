function [opdegs,maxdegs,mindegs] = getdeg(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [opdegs,maxdegs,mindegs] = getdeg(P) extract the degrees of the variables
% appearing in each parameter of a dopvar2d object P, with components 
%   R00 , R0x , R0y , R02;
%   Rx0 , Rxx , Rxy , Rx2;
%   Ry0 , Ryx , Ryy , Ry2;
%   R20 , R2x , R2y , R22;
% 
% INPUT
%   P: dopvar2d class object.
%
% OUTPUT
%   opdegs: struct with same fieldnames as P, with each field an nx4 array
%           containing the degrees of the variables in the associated
%           component in P. Each column corresponds to one of the four
%           variables, arranged as 
%               [P.var1(1), P.var2(1), P.var1(2), P.var2(2)];
%
%   maxdegs: struct with same fieldnames as P, with each field a 4x4 array
%           containing the maximal degrees of the variables in the
%           associated component in P. Maximal joint degrees in
%           combinations of variables are also provided, ordered as
%
%                 0 |         ss2 |         tt2 |         ss2*tt2
%          ---------|-------------|-------------|-----------------
%               ss1 |     ss1*ss2 |     ss1*tt2 |     ss1*ss2*tt2
%          ---------|-------------|-------------|-----------------
%               tt1 |     tt1*ss2 |     tt1*tt2 |     tt1*ss2*tt2 
%          ---------|-------------|-------------|-----------------
%           ss1*tt1 | ss1*tt1*ss2 | ss1*tt1*tt2 | ss1*tt1*ss2*tt2  
%
%           where ss1=P.var1(1), ss2=P.var1(2), tt1=P.var2(1), tt2=P.var2(2);
%
%   mindegs: similar to maxdegs, only providing minimal (total) degrees
%           instead of maximal degress in each variable
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - getdeg
%
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
% Initial coding DJ - 04_11_2022

if ~isa(P,'dopvar2d')
    error('Input must be a dopvar2d object.');
end

% ss1
if isa(P.var1(1),'polynomial')
    var11 = P.var1(1).varname;
elseif isa(P.var1,'char')
    var11 = P.var1(1);
end
% tt1
if isa(P.var2(1),'polynomial')
    var12 = P.var2(1).varname;
elseif isa(P.var2,'char')
    var12 = P.var2(1);
end
% ss2
if isa(P.var1(2),'polynomial')
    var21 = P.var1(2).varname;
elseif isa(P.var1,'char')
    var21 = P.var1(2);
end
% tt2
if isa(P.var2(2),'polynomial')
    var22 = P.var2(2).varname;
elseif isa(P.var2,'char')
    var22 = P.var2(2);
end

% Initialize the degree objects
opdegs = struct();
maxdegs = struct();
mindegs = struct();

opdegs.vars = reshape([P.var1,P.var2],4,1);
maxdegs.vars = reshape([P.var1,P.var2],4,1);
mindegs.vars = reshape([P.var1,P.var2],4,1);
        
fset = {'R00','R0x','R0y','R02','Rx0','Rxx','Rxy','Rx2',...
        'Ry0','Ryx','Ryy','Ry2','R20','R2x','R2y','R22'};   % fieldnames of the dopvar2d object
% For each parameter, extract the degmat
for f=1:length(fset)
    PR = P.(fset{f});
    if isa(PR,'double') || isempty(PR)
        % If the object is not polynomial, degrees are all zero
        opdegs.(fset{f}) = zeros(0,4);
        maxdegs.(fset{f}) = zeros(4,4);
        mindegs.(fset{f}) = zeros(4,4);
        
    elseif isa(PR,'polynomial') || isa(PR,'dpvar')
        % Find which variable in the polynomial corresponds to each
        % variable of the opvar
        idx1 = find(strcmp(PR.varname,var11));  % ss1
        idx2 = find(strcmp(PR.varname,var12));  % tt1
        idx3 = find(strcmp(PR.varname,var21));  % ss2
        idx4 = find(strcmp(PR.varname,var22));  % tt2
        
        PRdegs = PR.degmat;                 % Extract degtable in current vars
        nmons = size(PRdegs,1);
        opdegs.(fset{f}) = sparse(zeros(nmons,4));  % Initialize new degmat in all 4 vars
        maxdegs.(fset{f}) = zeros(4,4);             % Initialize max degs array
        mindegs.(fset{f}) = zeros(4,4);             % Initialize max degs array
        
        % % For each variable, add the degrees to the appropriate column of
        % % the degmat
        if ~isempty(idx1)   % Degrees in ss1
            opdegs.(fset{f})(:,1) = PRdegs(:,idx1);
        end
        if ~isempty(idx2)   % Degrees in tt1
            opdegs.(fset{f})(:,2) = PRdegs(:,idx2);
        end
        if ~isempty(idx3)   % Degrees in ss2
            opdegs.(fset{f})(:,3) = PRdegs(:,idx3);
        end
        if ~isempty(idx4)   % Degrees in tt2
            opdegs.(fset{f})(:,4) = PRdegs(:,idx4);
        end
        for l = 2:2^4   % Loop over all possible combinations of variables
            indcs = cell(1,4);
            [indcs{:}] = ind2sub(2*ones(1,4),l);        % Indices associated to position l in a 2x2x2x2 array
            log_indcs = cell2mat(indcs)==2;             % Logical values indicating which vars contribute to the considered joint degree
            maxdegs.(fset{f})(l) = max(sum(opdegs.(fset{f})(:,log_indcs),2)); % Take maximal joint degree in the considered variable
            mindegs.(fset{f})(l) = min(sum(opdegs.(fset{f})(:,log_indcs),2)); % Take maximal joint degree in the considered variable
        end
        
    elseif isa(PR,'cell')
        for j=1:numel(PR)
            PRR = PR{j};
            
            if isa(PRR,'double') || isempty(PRR)
                % If the object is not polynomial, degrees are all zero
                opdegs.(fset{f}){j} = zeros(0,4);
                maxdegs.(fset{f}){j} = zeros(4,4);
                mindegs.(fset{f}){j} = zeros(4,4);
            
            elseif isa(PRR,'polynomial') || isa(PRR,'dpvar')
                % Find which variable in the polynomial corresponds to each
                % variable of the opvar
                idx1 = find(strcmp(PRR.varname,var11));  % ss1
                idx2 = find(strcmp(PRR.varname,var12));  % tt1
                idx3 = find(strcmp(PRR.varname,var21));  % ss2
                idx4 = find(strcmp(PRR.varname,var22));  % tt2
                
                PRdegs = PRR.degmat;                 % Extract degtable in current vars
                nmons = size(PRdegs,1);
                opdegs.(fset{f}){j} = zeros(nmons,4); % Initialize new degmat in all 4 vars
                maxdegs.(fset{f}){j} = zeros(4,4);     % Initialize max degs array
                mindegs.(fset{f}){j} = zeros(4,4);     % Initialize max degs array
                
                % % For each variable, add the degrees to the appropriate column of
                % % the degmat
                if ~isempty(idx1)   % Degrees in ss1
                    opdegs.(fset{f}){j}(:,1) = PRdegs(:,idx1);
                end
                if ~isempty(idx2)   % Degrees in tt1
                    opdegs.(fset{f}){j}(:,2) = PRdegs(:,idx2);
                end
                if ~isempty(idx3)   % Degrees in ss2
                    opdegs.(fset{f}){j}(:,3) = PRdegs(:,idx3);
                end
                if ~isempty(idx4)   % Degrees in tt2
                    opdegs.(fset{f}){j}(:,4) = PRdegs(:,idx4);
                end
                for l = 2:2^4   % Loop over all possible combinations of variables
                    indcs = cell(1,4);
                    [indcs{:}] = ind2sub(2*ones(1,4),l);        % Indices associated to position l in a 2x2x2x2 array
                    log_indcs = cell2mat(indcs)==2;             % Logical values indicating which vars contribute to the considered joint degree
                    maxdegs.(fset{f}){j}(l) = max(sum(opdegs.(fset{f}){j}(:,log_indcs),2)); % Take maximal joint degree in the considered variable
                    mindegs.(fset{f}){j}(l) = min(sum(opdegs.(fset{f}){j}(:,log_indcs),2)); % Take maximal joint degree in the considered variable
                end
                
            else
                error(['Component ',fset{f},'{',num2str(j),'} of dopvar2d object is not properly specified'])
            end
        end
        
    else
        error(['Component ',fset{f},' of dopvar2d object is not properly specified'])
    end
end

opdegs.R22 = reshape(opdegs.R22,[3,3]);
maxdegs.R22 = reshape(maxdegs.R22,[3,3]);
mindegs.R22 = reshape(mindegs.R22,[3,3]);

end
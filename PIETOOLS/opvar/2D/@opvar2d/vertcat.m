function [Pcat] = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = vertcat(varargin) takes n inputs and concatentates them vertically,
% provided they satisfy the following criteria:
% 1) All inputs must of type 'opvar2d';
% 2) The input dimensions varargin{j}.dim(:,1) of all objects must match;
% 3) The spatial variables varargin{j}.var1 and varargin{j}.var2, as well 
%       as the domain varargin{j}.I of all objects must match;
% 4) Concatenation should make sense within the context of the opvar2d class,
%       that is, opvar2ds map to RxL2[x]xL2[y]xL2[x,y], 
%       so if the concatenated operator maps to any other space, 
%       e.g. L2[x]xRxL2[x,y]xL2[y], concatenation will be prohibited.
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - vertcat
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
%
% Initial coding DJ - 02_04_2021
%   ^ Based heavily on "@opvar"-vertcat code by SS ^
% DJ, 09/30/2023: Prohibit "ambiguous" concatenations;
% DJ, 12/15/2024: Allow conctanetation with matrix, polynomial, and dpvar;

% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% Extract the operators
a = varargin{1};    b = varargin{2};

% Check that inputs are of suitable type.
if isa(a,'dpvar')                                                           % DJ, 12/15/2024
    % Use dopvar concatenation.
    b = opvar2dopvar2d(b);
    Pcat = vertcat(vertcat(a,b),varargin{3:end});
    return
elseif isa(b,'dpvar')                                                       % DJ, 12/15/2024
    % Use dopvar concatenation.
    a = opvar2dopvar2d(a);
    Pcat = vertcat(vertcat(a,b),varargin{3:end});
    return
elseif isa(a,'double') || isa(a,'polynomial')                               % DJ, 12/15/2024
    if size(a,2)~=size(b,2)
        error('Cannot concatenate vertically: Input dimensions of operators do not match')
    end
    % Convert a to opvar2d, assuming to map to R, and with same input
    % dimension as b.
    a = mat2opvar(a,[[size(a,2);0;0;0],b.dim(:,2)],[b.var1,b.var2],b.I);
elseif isa(b,'double') || isa(b,'polynomial')                               % DJ, 12/15/2024
    if size(a,2)~=size(b,2)
        error('Cannot concatenate vertically: Input dimensions of operators do not match')
    end
    % Convert b to opvar2d, assuming to map to R, and with same input
    % dimension as a.
    b = mat2opvar(b,[[size(b,2);0;0;0],a.dim(:,2)],[a.var1,a.var2],a.I);
elseif ~isa(a,'opvar2d') || ~isa(b,'opvar2d')
    error('Concatenation of opvar2d and non-opvar2d objects is ambiguous, and currently not supported')
end
% Check that domain and variables match
if any(any(a.I~=b.I))|| ~all(strcmp(a.var1.varname,b.var1.varname)) || ~all(strcmp(a.var2.varname,b.var2.varname))
    error('Operators being concatenated have different spatial domains or variables');
end
% Check that the input dimensions match
a.dim = a.dim;  b.dim = b.dim;
if any(a.dim(:,2)~=b.dim(:,2))
    error('Cannot concatenate vertically: Input dimensions of opvar2d objects do not match')
end

% Avoid "ambiguous" concatenations
intype_a = find(a.dim(:,1),1,'last');   % What is the largest function space a maps to?
intype_b = find(b.dim(:,1),1,'first');  % What is the smallest function space b maps to?
if ~isempty(intype_a) && ~isempty(intype_b) && intype_b<intype_a
    % a has columns mapping from L2, but b has columns mapping from R.
    % Concatenation would place those columns of a to the right of
    % those columns of b in the opvar, which we currently prohibit...
    error('Proposed opvar2d concatenation is ambiguous, and currently prohibited')
end

% Finally, let's actually concatenate
% Initialize the concatenated operator
newdim = [a.dim(:,1)+b.dim(:,1),a.dim(:,2)];
Pcat = opvar2d([],newdim,a.I,a.var1,a.var2);

% Only concatenate columns which are nonempty
fset = {};
perform_cat = zeros(4,1);
if Pcat.dim(1,2)~=0
    fset = [fset,'R00','Rx0','Ry0','R20'];
    perform_cat(1) = 1;
end
if Pcat.dim(2,2)~=0
    fset = [fset,'R0x','Ryx'];
    perform_cat(2) = 1;
end
if Pcat.dim(3,2)~=0
    fset = [fset,'R0y','Rxy'];
    perform_cat(3) = 1;
end
if Pcat.dim(4,2)~=0
    fset = [fset,'R02'];
    perform_cat(4) = 1;
end

% Perform the concatenation
for f=fset
    Pcat.(f{:}) = [a.(f{:}); b.(f{:})];
end
for i=1:3
    if perform_cat(2)
        Pcat.Rxx{i,1} = [a.Rxx{i,1}; b.Rxx{i,1}];
        Pcat.Rx2{i,1} = [a.Rx2{i,1}; b.Rx2{i,1}];
    end
    if perform_cat(3)
        Pcat.Ryy{1,i} = [a.Ryy{1,i}; b.Ryy{1,i}];
        Pcat.Ry2{1,i} = [a.Ry2{1,i}; b.Ry2{1,i}];
    end
    if perform_cat(4)
        Pcat.R2x{i,1} = [a.R2x{i,1}; b.R2x{i,1}];
        Pcat.R2y{1,i} = [a.R2y{1,i}; b.R2y{1,i}];
        for j=1:3
            Pcat.R22{i,j} = [a.R22{i,j}; b.R22{i,j}];
        end
    end
end

% For concatenation of more than two objects, just repeat
if nargin>2 
    Pcat = vertcat(Pcat, varargin{3:end});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Old version, that allows "ambiguous" concatenation   %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = vertcat(varargin) takes n-inputs and concatentates them vertically,
% provided they satisfy the following criterias.
% 1) Atleast one input is an opvar2d variable.
% 2) If all the inputs are not opvar2d, then the operator maps from R to
% RxL2 or L2 to L2. 
% 3) Currently, it supports RxL2 to RxL2 concatenation only if ALL the inputs are
% opvar2d.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% if nargin==1
%     Pcat = varargin{1};
% else
%     a = varargin{1};
%     b = varargin{2};
%     
%     if isa(a,'opvar2d') % components must have consistent dimensions
%         a.dim = a.dim;
%     end
%     if isa(b,'opvar2d')
%         b.dim = b.dim;
%     end
%     
%     if ~isa(a,'opvar2d') 
%         if ~isa(b,'opvar2d')
%             if size(a,2)~=size(b,2)
%                 error('Cannot concatentate vertically. A and B have different input dimensions');
%             end
%             Pcat = [a;b];
%         else
%             bdim = b.dim;                
%             if size(a,2)~=sum(bdim(:,2))
%                 error("Cannot concatentate vertically. A and B have different input dimensions");
%             end
%             Pcat = b;
%             Pcat.dim = [b.dim(:,1)+(b.dim(:,1)~=0).*size(a,1),b.dim(:,2)];
%             if ~isempty(a)
%                 if all(bdim(2:4,2) == 0) % a() is from R to R
%                     Pcat.R00 = [a; b.R00]; 
%                     Pcat.R0x = [zeros(size(a,1),b.dim(2,2)); b.R0x];
%                     Pcat.R0y = [zeros(size(a,1),b.dim(3,2)); b.R0y];
%                     Pcat.R02 = [zeros(size(a,1),b.dim(4,2)); b.R02];
%                 elseif all(bdim([1,3:4],2) == 0) % a() is from L2[s1] to L2[s1], Note: For L2 to R, a must be opvar2d 
%                     Pcat.Rx0 = [zeros(size(a,1),b.dim(1,2)); b.Rx0];
%                     Pcat.Rxx{1,1} = [a; b.Rxx{1,1}];
%                     Pcat.Rxx{2,1} = [zeros(size(a,1),b.dim(2,2)); b.Rxx{2,1}];
%                     Pcat.Rxx{3,1} = [zeros(size(a,1),b.dim(2,2)); b.Rxx{3,1}];
%                     Pcat.Rxy = [zeros(size(a,1),b.dim(3,2)); b.Rxy];
%                     Pcat.Rx2{1,1} = [zeros(size(a,1),b.dim(4,2)); b.Rx2{1,1}];
%                     Pcat.Rx2{2,1} = [zeros(size(a,1),b.dim(4,2)); b.Rx2{2,1}];
%                     Pcat.Rx2{3,1} = [zeros(size(a,1),b.dim(4,2)); b.Rx2{3,1}];
%                 elseif all(bdim([1:2,4],2) == 0) % a() is from L2[s2] to L2[s2], Note: For L2 to R, a must be opvar2d 
%                     Pcat.Ry0 = [zeros(size(a,1),b.dim(1,2)); b.Ry0];
%                     Pcat.Ryx = [zeros(size(a,1),b.dim(2,2)); b.Ryx];
%                     Pcat.Ryy{1,1} = [a; b.Ryy{1,1}];
%                     Pcat.Ryy{1,2} = [zeros(size(a,1),b.dim(3,2)); b.Ryy{1,2}];
%                     Pcat.Ryy{1,3} = [zeros(size(a,1),b.dim(3,2)); b.Ryy{1,3}];
%                     Pcat.Ry2{1,1} = [zeros(size(a,1),b.dim(4,2)); b.Ry2{1,1}];
%                     Pcat.Ry2{1,2} = [zeros(size(a,1),b.dim(4,2)); b.Ry2{1,2}];
%                     Pcat.Ry2{1,3} = [zeros(size(a,1),b.dim(4,2)); b.Ry2{1,3}];
%                 elseif all(bdim(1:3,2) == 0) % a() is from L2[s1,s2] to L2[s1,s2], Note: For L2 to R, a must be opvar2d 
%                     Pcat.R20 = [zeros(size(a,1),b.dim(1,2)); b.R20];
%                     Pcat.R2x{1,1} = [zeros(size(a,1),b.dim(2,2)); b.R2x{1,1}];
%                     Pcat.R2x{2,1} = [zeros(size(a,1),b.dim(2,2)); b.R2x{2,1}];
%                     Pcat.R2x{3,1} = [zeros(size(a,1),b.dim(2,2)); b.R2x{3,1}];
%                     Pcat.R2y{1,1} = [zeros(size(a,1),b.dim(3,2)); b.R2y{1,1}];
%                     Pcat.R2y{1,2} = [zeros(size(a,1),b.dim(3,2)); b.R2y{1,2}];
%                     Pcat.R2y{1,3} = [zeros(size(a,1),b.dim(3,2)); b.R2y{1,3}];
%                     Pcat.R22{1,1} = [a; b.R22{1,1}];
%                     Pcat.R22{2,1} = [zeros(size(a,1),b.dim(4,2)); b.R22{2,1}];
%                     Pcat.R22{3,1} = [zeros(size(a,1),b.dim(4,2)); b.R22{3,1}];
%                     Pcat.R22{1,2} = [zeros(size(a,1),b.dim(4,2)); b.R22{1,2}];
%                     Pcat.R22{2,2} = [zeros(size(a,1),b.dim(4,2)); b.R22{2,2}];
%                     Pcat.R22{3,2} = [zeros(size(a,1),b.dim(4,2)); b.R22{3,2}];
%                     Pcat.R22{1,3} = [zeros(size(a,1),b.dim(4,2)); b.R22{1,3}];
%                     Pcat.R22{2,3} = [zeros(size(a,1),b.dim(4,2)); b.R22{2,3}];
%                     Pcat.R22{3,3} = [zeros(size(a,1),b.dim(4,2)); b.R22{3,3}];
%                 else %find if such an operation is valid in any useful scenario and implement it
%                     error('Cannot concatenate vertically. This feature is not yet supported.');
%                 end
%             end
%         end
%     elseif ~isa(b,'opvar2d') 
%         adim = a.dim;
%         if size(b,2)~=sum(adim(:,2))
%             error("Cannot concatentate vertically. A and B have different input dimensions");
%         end
%         Pcat = a;
%         Pcat.dim = [a.dim(:,1)+(a.dim(:,1)~=0).*size(b,1),a.dim(:,2)];
%         if ~isempty(b)
%             if all(adim(2:4,2) == 0) % b() is from R to R
%                 Pcat.R00 = [a.R00; b]; 
%                 Pcat.R0x = [a.R0x; zeros(size(b,1),a.dim(2,2))];
%                 Pcat.R0y = [a.R0y; zeros(size(b,1),a.dim(3,2))];
%                 Pcat.R02 = [a.R02; zeros(size(b,1),a.dim(4,2))];
%             elseif all(adim([1,3:4],2) == 0) % b() is from L2[s1] to L2[s1]
%                 Pcat.Rx0 = [a.Rx0; zeros(size(b,1),a.dim(1,2))];
%                 Pcat.Rxx{1,1} = [a.Rxx{1,1}; b];
%                 Pcat.Rxx{2,1} = [a.Rxx{2,1}; zeros(size(b,1),a.dim(2,2))];
%                 Pcat.Rxx{3,1} = [a.Rxx{3,1}; zeros(size(b,1),a.dim(2,2))];
%                 Pcat.Rxy = [a.Rxy; zeros(size(b,1),a.dim(3,2))];
%                 Pcat.Rx2{1,1} = [a.Rx2{1,1}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.Rx2{2,1} = [a.Rx2{2,1}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.Rx2{3,1} = [a.Rx2{3,1}; zeros(size(b,1),a.dim(4,2))];
%             elseif all(adim([1:2,4],2) == 0) % b() is from L2[s2] to L2[s2]
%                 Pcat.Ry0 = [a.Ry0; zeros(size(b,1),a.dim(1,2))];
%                 Pcat.Ryx = [a.Ryx; zeros(size(b,1),a.dim(2,2))];
%                 Pcat.Ryy{1,1} = [a.Ryy{1,1}; b];
%                 Pcat.Ryy{1,2} = [a.Ryy{1,2}; zeros(size(b,1),a.dim(3,2))];
%                 Pcat.Ryy{1,3} = [a.Ryy{1,3}; zeros(size(b,1),a.dim(3,2))];
%                 Pcat.Ry2{1,1} = [a.Ry2{1,1}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.Ry2{1,2} = [a.Ry2{1,2}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.Ry2{1,3} = [a.Ry2{1,3}; zeros(size(b,1),a.dim(4,2))];
%             elseif all(adim(1:3,2) == 0) % b() is from L2[s1,s2] to L2[s1,s2]
%                 Pcat.R20 = [a.R20; zeros(size(b,1),a.dim(1,2))];
%                 Pcat.R2x{1,1} = [a.R2x{1,1}; zeros(size(b,1),a.dim(2,2))];
%                 Pcat.R2x{2,1} = [a.R2x{2,1}; zeros(size(b,1),a.dim(2,2))];
%                 Pcat.R2x{3,1} = [a.R2x{3,1}; zeros(size(b,1),a.dim(2,2))];
%                 Pcat.R2y{1,1} = [a.R2y{1,1}; zeros(size(b,1),a.dim(3,2))];
%                 Pcat.R2y{1,2} = [a.R2y{1,2}; zeros(size(b,1),a.dim(3,2))];
%                 Pcat.R2y{1,3} = [a.R2y{1,3}; zeros(size(b,1),a.dim(3,2))];
%                 Pcat.R22{1,1} = [a.R22{1,1}; b];
%                 Pcat.R22{2,1} = [a.R22{2,1}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.R22{3,1} = [a.R22{3,1}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.R22{1,2} = [a.R22{1,2}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.R22{2,2} = [a.R22{2,2}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.R22{3,2} = [a.R22{3,2}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.R22{1,3} = [a.R22{1,3}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.R22{2,3} = [a.R22{2,3}; zeros(size(b,1),a.dim(4,2))];
%                 Pcat.R22{3,3} = [a.R22{3,3}; zeros(size(b,1),a.dim(4,2))];
%             else %find if such an operation is valid in any useful scenario and implement it
%                 error('Cannot concatenate vertically. This feature is not yet supported.');
%             end
%         end
%     else
%         if any(b.dim(:,2)~=a.dim(:,2))
%             error('Cannot concatentate horizontally. A and B have different column dimensions');
%         elseif any(any(b.I~=a.I))
%             error('Cannot concatentate horizontally: A and B have different domains');
%         end
%         % Initialize the concatenated operator
%         newdim = [a.dim(:,1)+b.dim(:,1),a.dim(:,2)];
%         Pcat = opvar2d([],newdim,a.I,a.var1,a.var2);
%         
%         % Only concatenate columns which are nonempty
%         fset = {};
%         c = zeros(4,1);
%         if Pcat.dim(1,2)~=0
%             fset = [fset,'R00','Rx0','Ry0','R20'];
%             c(1) = 1;
%         end
%         if Pcat.dim(2,2)~=0
%             fset = [fset,'R0x','Ryx'];
%             c(2) = 1;
%         end
%         if Pcat.dim(3,2)~=0
%             fset = [fset,'R0y','Rxy'];
%             c(3) = 1;
%         end
%         if Pcat.dim(4,2)~=0
%             fset = [fset,'R02'];
%             c(4) = 1;
%         end
%         
%         % Perform the concatenation
%         for f=fset
%             Pcat.(f{:}) = [a.(f{:}); b.(f{:})];
%         end
%         for i=1:3
%             if c(2)
%                 Pcat.Rxx{i,1} = [a.Rxx{i,1}; b.Rxx{i,1}];
%                 Pcat.R2x{i,1} = [a.R2x{i,1}; b.R2x{i,1}];
%             end
%             if c(3)
%                 Pcat.Ryy{1,i} = [a.Ryy{1,i}; b.Ryy{1,i}];
%                 Pcat.R2y{1,i} = [a.R2y{1,i}; b.R2y{1,i}];
%             end
%             if c(4)
%                 Pcat.Rx2{i,1} = [a.Rx2{i,1}; b.Rx2{i,1}];
%                 Pcat.Ry2{1,i} = [a.Ry2{1,i}; b.Ry2{1,i}];
%                 for j=1:3
%                     Pcat.R22{i,j} = [a.R22{i,j}; b.R22{i,j}];
%                 end
%             end
%         end
%     end
%     if nargin>2 % Continue concatenation if inputs are more than 2
%         Pcat = vertcat(Pcat, varargin{3:end});
%     end
% end
% end
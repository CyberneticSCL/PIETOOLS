function [Pcat] = horzcat(varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = horzcat(varargin) takes n inputs and concatentates them horizontally,
% provided they satisfy the following criteria:
% 1) All inputs must of type 'opvar2d';
% 2) The output dimensions varargin{j}.dim(:,2) of all objects must match;
% 3) The spatial variables varargin{j}.var1 and varargin{j}.var2, as well 
%       as the domain varargin{j}.I of all objects must match;
% 4) Concatenation should make sense within the context of the opvar2d class,
%       that is, opvar2ds take inputs from RxL2[x]xL2[y]xL2[x,y], 
%       so if the concatenated operator maps any other space, 
%       e.g. L2[x]xRxL2[x,y]xL2[y], concatenation will be prohibited.
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - horzcat
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
%   ^ Based heavily on "@opvar"-horzcat code by SS ^
% DJ - 09/30/23: Prohibit "ambiguous" concatenations.

% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% Extract the operators
a = varargin{1};    b = varargin{2};

% Currently support only opvar-opvar concatenation
if ~isa(a,'opvar2d') || ~isa(b,'opvar2d')
    error('Concatenation of opvar2d and non-opvar2d objects is ambiguous, and currently not supported')
end
% Check that domain and variables match
if any(any(a.I~=b.I))|| ~all(strcmp(a.var1.varname,b.var1.varname)) || ~all(strcmp(a.var2.varname,b.var2.varname))
    error('Operators being concatenated have different spatial domains or variables');
end
% Check that the output dimensions match
a.dim = a.dim;  b.dim = b.dim;
if any(a.dim(:,1)~=b.dim(:,1))
    error('Cannot concatenate horizontally: Output dimensions of opvar2d objects do not match')
end

% Avoid "ambiguous" concatenations
intype_a = find(a.dim(:,2),1,'last');   % What is the largest function space a maps?
intype_b = find(b.dim(:,2),1,'first');  % What is the smallest function space b maps?
if ~isempty(intype_a) && ~isempty(intype_b) && intype_b<intype_a
    % a has columns mapping from L2, but b has columns mapping from R.
    % Concatenation would place those columns of a to the right of
    % those columns of b in the opvar, which we currently prohibit...
    error('Proposed opvar2d concatenation is ambiguous, and currently prohibited')
end

% Finally, let's actually concatenate
% Initialize the concatenated operator
newdim = [a.dim(:,1),a.dim(:,2)+b.dim(:,2)];
Pcat = opvar2d([],newdim,a.I,a.var1,a.var2);

% Only concatenate rows which are nonempty
fset = {};
perform_cat = zeros(4,1);
if Pcat.dim(1,1)~=0
    fset = [fset,'R00','R0x','R0y','R02'];
    perform_cat(1) = 1;
end
if Pcat.dim(2,1)~=0
    fset = [fset,'Rx0','Rxy'];
    perform_cat(2) = 1;
end
if Pcat.dim(3,1)~=0
    fset = [fset,'Ry0','Ryx'];
    perform_cat(3) = 1;
end
if Pcat.dim(4,1)~=0
    fset = [fset,'R20'];
    perform_cat(4) = 1;
end

% Perform the concatenation
for f=fset
    Pcat.(f{:}) = [a.(f{:}) b.(f{:})];
end
for i=1:3
    if perform_cat(2)
        Pcat.Rxx{i,1} = [a.Rxx{i,1} b.Rxx{i,1}];
        Pcat.Rx2{i,1} = [a.Rx2{i,1} b.Rx2{i,1}];
    end
    if perform_cat(3)
        Pcat.Ryy{1,i} = [a.Ryy{1,i} b.Ryy{1,i}];
        Pcat.Ry2{1,i} = [a.Ry2{1,i} b.Ry2{1,i}];
    end
    if perform_cat(4)
        Pcat.R2x{i,1} = [a.R2x{i,1} b.R2x{i,1}];
        Pcat.R2y{1,i} = [a.R2y{1,i} b.R2y{1,i}];
        for j=1:3
            Pcat.R22{i,j} = [a.R22{i,j} b.R22{i,j}];
        end
    end
end

% For concatenation of more than two objects, just repeat
if nargin>2 
    Pcat = horzcat(Pcat, varargin{3:end});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Old version, that allows "ambiguous" concatenation   %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = horzcat(varargin) takes n-inputs and concatentates them horizontally,
% provided they satisfy the following criterias.
% 1) Atleast one input is an opvar2d variable.
% 2) If all the inputs are not opvar2d, then the operator maps from RxL2 to
%    L2 or R to R. 
% 3) Currently, it supports RxL2 to RxL2 concatenation only if ALL the inputs are
%    opvar2d.
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
%     if isa(a,'opvar2d')
%         a.dim = a.dim;
%     end
%     if isa(b,'opvar2d')
%         b.dim = b.dim;
%     end
%     
%     if ~isa(a,'opvar2d')
%         if ~isa(b,'opvar2d') % both are not opvar2d variables
%             if size(a,1)~=size(b,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             Pcat = [a b];
%         elseif all(b.dim(2:4,1)==0) % a() is from R to R, Note: For L2 to R, a needs to be an opvar2d
%             if size(a,1)~=b.dim(1,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end 
%             Pcat = b;
%             Pcat.dim(1,2) = Pcat.dim(1,2) + size(a,2);
%             if ~isempty(a)
%                 Pcat.R00 = [a b.R00];
%                 Pcat.Rx0 = [zeros(b.dim(2,1),size(a,2)) b.Rx0];
%                 Pcat.Ry0 = [zeros(b.dim(3,1),size(a,2)) b.Ry0];
%                 Pcat.R20 = [zeros(b.dim(4,1),size(a,2)) b.R20];
%             end
%         elseif all(b.dim([1,3:4],1)==0) % a() is from R to L2[s1], Note: For L2 to L2, a needs to be an opvar2d
%             if size(a,1)~=b.dim(2,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end 
%             Pcat = b;
%             Pcat.dim(2,2) = Pcat.dim(2,2) + size(a,2);
%             if ~isempty(a)
%                 Pcat.R00 = [zeros(b.dim(1,1),size(a,2)) b.R00];
%                 Pcat.Rx0 = [a b.Rx0];
%                 Pcat.Ry0 = [zeros(b.dim(3,1),size(a,2)) b.Ry0];
%                 Pcat.R20 = [zeros(b.dim(4,1),size(a,2)) b.R20];
%             end
%         elseif all(b.dim([1:2,4],1)==0) % a() is from R to L2[s2], Note: For L2 to L2, a needs to be an opvar2d
%             if size(a,1)~=b.dim(3,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end 
%             Pcat = b;
%             Pcat.dim(3,2) = Pcat.dim(3,2) + size(a,2);
%             if ~isempty(a)
%                 Pcat.R00 = [zeros(b.dim(1,1),size(a,2)) b.R00];
%                 Pcat.Rx0 = [zeros(b.dim(2,1),size(a,2)) b.Ry0];
%                 Pcat.Ry0 = [a b.Ry0];
%                 Pcat.R20 = [zeros(b.dim(4,1),size(a,2)) b.R20];
%             end
%         elseif all(b.dim(1:3,1)==0) % a() is from R to L2[s1,s2], Note: For L2 to L2, a needs to be an opvar2d
%             if size(a,1)~=b.dim(4,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end 
%             Pcat = b;
%             Pcat.dim(4,2) = Pcat.dim(4,2) + size(a,2);
%             if ~isempty(a)
%                 Pcat.R00 = [zeros(b.dim(1,1),size(a,2)) b.R00];
%                 Pcat.Rx0 = [zeros(b.dim(2,1),size(a,2)) b.Rx0];
%                 Pcat.Ry0 = [zeros(b.dim(3,1),size(a,2)) b.Ry0];
%                 Pcat.R20 = [a b.R20];
%             end
%         else %find if such an operation is valid in any useful scenario and implement it
%             error('Cannot concatenate horizontally. This feature is not yet supported.');
%         end
%     elseif ~isa(b,'opvar2d')
%         if all(a.dim(2:4,1)==0) % a() is from R to R, Note: L2 to L2 needs to be an opvar
%             if size(b,1)~=a.dim(1,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             Pcat = a;
%             Pcat.dim(1,2) = Pcat.dim(1,2) + size(b,2);
%             if ~isempty(b)
%                 Pcat.R00 = [a.R00 b];
%                 Pcat.Rx0 = [a.Rx0 zeros(a.dim(2,1),size(b,2))];
%                 Pcat.Ry0 = [a.Ry0 zeros(a.dim(3,1),size(b,2))];
%                 Pcat.R20 = [a.R20 zeros(a.dim(4,1),size(b,2))];
%             end
%         elseif all(a.dim([1,3:4],1)==0) % a() is from R to L2[s1], Note: L2 to L2 needs to be an opvar
%             if size(b,1)~=a.dim(2,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             Pcat = a;
%             Pcat.dim(2,2) = Pcat.dim(2,2) + size(b,2);
%             if ~isempty(b)
%                 Pcat.R00 = [a.R00 zeros(a.dim(1,1),size(b,2))];
%                 Pcat.Rx0 = [a.Rx0 b];
%                 Pcat.Ry0 = [a.Ry0 zeros(a.dim(3,1),size(b,2))];
%                 Pcat.R20 = [a.R20 zeros(a.dim(4,1),size(b,2))];
%             end
%         elseif all(a.dim([1:2,4],1)==0) % a() is from R to L2[s2], Note: L2 to L2 needs to be an opvar
%             if size(b,1)~=a.dim(3,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             Pcat = a;
%             Pcat.dim(3,2) = Pcat.dim(3,2) + size(b,2);
%             if ~isempty(b)
%                 Pcat.R00 = [a.R00 zeros(a.dim(1,1),size(b,2))];
%                 Pcat.Rx0 = [a.Rx0 zeros(a.dim(2,1),size(b,2))];
%                 Pcat.Ry0 = [a.Ry0 b];
%                 Pcat.R20 = [a.R20 zeros(a.dim(4,1),size(b,2))];
%             end
%         elseif all(a.dim(1:3,1)==0) % a() is from R to L2[s1,s2], Note: L2 to L2 needs to be an opvar
%             if size(b,1)~=a.dim(4,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             Pcat = a;
%             Pcat.dim(4,2) = Pcat.dim(4,2) + size(b,2);
%             if ~isempty(b)
%                 Pcat.R00 = [a.R00 zeros(a.dim(1,1),size(b,2))];
%                 Pcat.Rx0 = [a.Rx0 zeros(a.dim(2,1),size(b,2))];
%                 Pcat.Ry0 = [a.Ry0 zeros(a.dim(3,1),size(b,2))];
%                 Pcat.R20 = [a.R20 b];
%             end
%         else %find if such a operation is valid in any useful scenario and implement it
%             error('Cannot concatenate horizontally. This feature is not yet supported.');
%         end
%     else % Both arguments are opvar2d
%         if any(b.dim(:,1)~=a.dim(:,1))
%             error('Cannot concatentate horizontally. A and B have different row dimensions');
%         elseif any(any(b.I~=a.I))
%             error('Cannot concatentate horizontally: A and B have different domains');
%         end
%         % Initialize the concatenated operator
%         newdim = [a.dim(:,1),a.dim(:,2)+b.dim(:,2)];
%         Pcat = opvar2d([],newdim,a.I,a.var1,a.var2);
%         
%         % Only concatenate rows which are nonempty
%         fset = {};
%         r = zeros(4,1);
%         if Pcat.dim(1,1)~=0
%             fset = [fset,'R00','R0x','R0y','R02'];
%             r(1) = 1;
%         end
%         if Pcat.dim(2,1)~=0
%             fset = [fset,'Rx0','Rxy'];
%             r(2) = 1;
%         end
%         if Pcat.dim(3,1)~=0
%             fset = [fset,'Ry0','Ryx'];
%             r(3) = 1;
%         end
%         if Pcat.dim(4,1)~=0
%             fset = [fset,'R20'];
%             r(4) = 1;
%         end
%         
%         % Perform the concatenation
%         for f=fset
%             Pcat.(f{:}) = [a.(f{:}) b.(f{:})];
%         end
%         for i=1:3
%             if r(2)
%                 Pcat.Rxx{i,1} = [a.Rxx{i,1} b.Rxx{i,1}];
%                 Pcat.Rx2{i,1} = [a.Rx2{i,1} b.Rx2{i,1}];
%             end
%             if r(3)
%                 Pcat.Ryy{1,i} = [a.Ryy{1,i} b.Ryy{1,i}];
%                 Pcat.Ry2{1,i} = [a.Ry2{1,i} b.Ry2{1,i}];
%             end
%             if r(4)
%                 Pcat.R2x{i,1} = [a.R2x{i,1} b.R2x{i,1}];
%                 Pcat.R2y{1,i} = [a.R2y{1,i} b.R2y{1,i}];
%                 for j=1:3
%                     Pcat.R22{i,j} = [a.R22{i,j} b.R22{i,j}];
%                 end
%             end
%         end
%     end
%     if nargin>2 % check if there are more than 2 objects that need to stacked
%         Pcat = horzcat(Pcat, varargin{3:end});
%     end
% end
% end
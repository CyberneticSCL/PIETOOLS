function [Pcat] = horzcat(varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = horzcat(varargin) takes n inputs and concatentates them horizontally,
% provided they satisfy the following criteria:
% 1) All inputs must of type 'opvar';
% 2) The output dimensions varargin{j}.dim(:,1) of all opvar objects must
%       match;
% 3) The spatial variables varargin{j}.var1 and varargin{j}.var2, as well 
%       as the domain varargin{j}.I of all opvar objects must match;
% 4) Concatenation should make sense within the context of the opvar class,
%       that is, opvars take inputs from RxL2, never from L2xR. We cannot
%       concatenate e.g. [A,B] for A:L2-->R and B:R-->R.
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - horzcat
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 7_26_2019
% SS - 8/1/19 added all cases of concatenation
% SS - 8/27/19 added dimension correction to work with new opvar class that
% has automatic dimension calculation
% Adjusted so that dpvar with opvar returns dopvar, DJ 12/30/2021.
% DJ - 09/30/23: Prohibit "ambiguous" concatenations.

% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% Extract the operators
a = varargin{1};    b = varargin{2};

% Currently support only opvar-opvar concatenation
if ~isa(a,'opvar') || ~isa(b,'opvar')
    error('Concatenation of opvar and non-opvar objects is ambiguous, and currently not supported')
end
% Check that domain and variables match
if any(a.I~=b.I)|| ~strcmp(a.var1.varname,b.var1.varname) || ~strcmp(a.var2.varname,b.var2.varname)
    error('Operators being concatenated have different intervals or different independent variables');
end
% Check that the output dimensions match
a.dim = a.dim;  b.dim = b.dim;  % correction to make components have consistent dimensions 8/27-ss
if any(a.dim(:,1)~=b.dim(:,1))
    error('Cannot concatenate horizontally: Output dimensions of opvar objects do not match')
end

% Avoid "ambiguous" concatenations
if (a.dim(2,2)~=0 && b.dim(1,2)~=0)
    % a has columns mapping from L2, but b has columns mapping from R.
    % Concatenation would place those columns of a to the right of
    % those columns of b in the opvar, which we currently prohibit...
    error('Proposed opvar concatenation is ambiguous, and currently prohibited')
end

% Finally, let's actually concatenate
opvar Pcat;
Pcat.I = a.I;   Pcat.var1 = a.var1;     Pcat.var2 = a.var2;
fset = {'P', 'Q1', 'Q2'};
for i=fset
    Pcat.(i{:}) = [a.(i{:}) b.(i{:})];
end
fset = {'R0','R1','R2'};
for i=fset
    Pcat.R.(i{:}) = [a.R.(i{:}) b.R.(i{:})];
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
% provided they satisfy the following criteria:
% 1) Atleast one input is an opvar variable.
% 2) If all the inputs are not opvar, then the operator maps from RxL2 to
%    L2 or R to R. 
% 3) Currently, it supports RxL2 to RxL2 concatenation only if ALL the inputs are
%    opvar.
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if nargin==1
%     Pcat = varargin{1};
% else
%     a = varargin{1};
%     b = varargin{2};
%     
%     opvar Pcat;
%     if isa(a,'opvar')
%         Pcat.I = a.I; Pcat.var1 = a.var1; Pcat.var2 = a.var2;
%     elseif isa(b,'opvar')
%         Pcat.I = b.I; Pcat.var1 = b.var1; Pcat.var2 = b.var2;
%     elseif isa(a,'opvar')&&isa(b,'opvar')
%         if any(a.I~=b.I)||(a.var1~=b.var1)||(a.var2~=b.var2)
%             error('Operators being concatenated have different intervals or different independent variables');
%         end
%     end
%     if isa(a,'opvar') % correction to make components have consistent dimensions 8/27-ss
%         a.dim = a.dim;
%     end
%     if isa(b,'opvar') % correction to make components have consistent dimensions 8/27-ss
%         b.dim = b.dim;
%     end
%     
%     if isa(a,'dpvar')  % DJ, 12/30/2021
%         b = opvar2dopvar(b);
%         Pcat = horzcat(a,b);
%     elseif isa(b,'dpvar')
%         a = opvar2dopvar(a);
%         Pcat = horzcat(a,b);
%     elseif ~isa(a,'opvar')
%         if ~isa(b,'opvar') % both are not opvar variables
%             if size(a,1)~=size(b,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             Pcat = [a b];
%         elseif b.dim(1,1)==0
%             % b:RxL2-->L2, assume a:R-->L2
%             % NOTE: a:L2-->L2 must be specified as opvar
%             if size(a,1)~=b.dim(2,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end 
%             %Pcat.P = [zeros(b.dim(1,1),size(a,2)), b.P];
%             Pcat.Q2 = [a,b.Q2];
%             Pcat.R.R0 = b.R.R0;
%             Pcat.R.R1 = b.R.R1;
%             Pcat.R.R2 = b.R.R2;
% %             Pcat.Q2 = b.Q2;
% %             Pcat.R.R0 = [a, b.R.R0];
% %             Pcat.R.R1 = [zeros(b.dim(2,1),size(a,2)), b.R.R1];
% %             Pcat.R.R2 = [zeros(b.dim(2,1),size(a,2)), b.R.R2];
%         elseif b.dim(2,1)==0 
%             % b:RxL2-->R, assume a:R-->R
%             if size(a,1)~=b.dim(1,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             if isa(a,'polynomial') || isa(b,'dpvar')
%                 try a=double(a);
%                 catch
%                     error('Convert arguments to opvar for this concatenation')
%                 end
%             end
%             Pcat.P = [a b.P];
%             Pcat.Q1 = b.Q1;
%             %Pcat.Q2 = [zeros(b.dim(2,1),size(a,2)) b.Q2];
%         else %find if such an operation is valid is any useful scenario and implement it
%             error('Cannot concatenate horizontally. This feature is not yet supported.');
%         end
%     elseif ~isa(b,'opvar')
%         if a.dim(1,1)==0
%             % a:RxL2-->L2, assume b:R-->L2
%             % NOTE: b:L2-->L2 must be specified as opvar
%             if size(b,1)~=a.dim(2,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             %Pcat.P = [a.P zeros(a.dim(1,1),size(b,2))];
%             Pcat.Q2 = [a.Q2, b];
%             Pcat.R.R0 = a.R.R0;
%             Pcat.R.R1 = a.R.R1;
%             Pcat.R.R2 = a.R.R2;
% %             Pcat.Q2 = a.Q2;
% %             Pcat.R.R0 = [a.R.R0, b];
% %             Pcat.R.R1 = [a.R.R1, zeros(a.dim(2,1),size(b,2))];
% %             Pcat.R.R2 = [a.R.R2, zeros(a.dim(2,1),size(b,2))];
%         elseif a.dim(2,1)==0
%             % a:RxL2-->R, assume b:R-->R
%             if size(b,1)~=a.dim(1,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
%             if isa(a,'polynomial') || isa(b,'dpvar')
%                 try a=double(a);
%                 catch
%                     error('Convert arguments to opvar for this concatenation')
%                 end
%             end
%             Pcat.P = [a.P b];
%             Pcat.Q1 = a.Q1;
%             %Pcat.Q2 = [a.Q2 zeros(a.dim(2,1),size(b,2))];
%         else %find if such a operation is valid is any useful scenario and implement it
%             error('Cannot concatenate horizontally.This feature is not yet supported.');
%         end
%     else
%         if any(b.dim(:,1)~=a.dim(:,1))
%             error('Cannot concatentate horizontally. A and B have different output dimensions');
%         end
%         fset = {'P', 'Q1', 'Q2'};
%         for i=fset
%             Pcat.(i{:}) = [a.(i{:}) b.(i{:})];
%         end
%         fset = {'R0','R1','R2'};
%         for i=fset
%             Pcat.R.(i{:}) = [a.R.(i{:}) b.R.(i{:})];
%         end
%     end
%     if nargin>2 % check if there are more than 2 objects that need to stacked
%         Pcat = horzcat(Pcat, varargin{3:end});
%     end
% end
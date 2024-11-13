function [Pcat] = horzcat(varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = horzcat(varargin) takes n inputs and concatentates them horizontally,
% provided they satisfy the following criteria:
% 1) At least one input must be of type 'dopvar', and all others must be of
%       type 'opvar' or 'dopvar';
% 2) The output dimensions varargin{j}.dim(:,1) of all objects must match;
% 3) The spatial variables varargin{j}.var1 and varargin{j}.var2, as well 
%       as the domain varargin{j}.I of all objects must match;
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
% SS - 9/26/2021 adding bug fixes to when concatenation has elements that
% are polynomials or just opvars
% DJ - 12/30/2021 Adjusted to assure opvar with dopvar returns dopvar
% DJ - 09/30/23: Prohibit "ambiguous" concatenations.
% DB - 11/13/24: "ambiguous" concatenations allowed, warnning added.

% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% Extract the operators
a = varargin{1};    b = varargin{2};

% Currently support only opvar-opvar concatenation
if (~isa(a,'opvar') && ~isa(a,'dopvar')) || (~isa(b,'opvar') && ~isa(b,'dopvar'))
    error('dopvar  objects can only be concatenated with opvar and dopvar objects')
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
    warnning('Proposed opvar concatenation is ambiguous')
end

% Finally, let's actually concatenate
dopvar Pcat;
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
% 1) Atleast one input is an dopvar variable.
% 2) If all the inputs are not dopvar, then the operator maps from RxL2 to
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
%     if isa(a,'dopvar') % correction to make components have consistent dimensions 8/27-ss
%         a.dim = a.dim;
%     end
%     if isa(b,'dopvar') % correction to make components have consistent dimensions 8/27-ss
%         b.dim = b.dim;
%     end
%     
%     dopvar Pcat;
%     if isa(a,'dopvar')
%         Pcat.I = a.I; Pcat.var1 = a.var1; Pcat.var2 = a.var2;
%     elseif isa(b,'dopvar')
%         Pcat.I = b.I; Pcat.var1 = b.var1; Pcat.var2 = b.var2;
%     elseif isa(a,'dopvar')&&isa(b,'dopvar')
%         if any(a.I~=b.I)||(a.var1~=b.var1)||(a.var2~=b.var2)
%             error('Operators being concatenated have different intervals or different independent variables');
%         end
%     end
%     
%     if ~isa(a,'dopvar')
%         if ~isa(b,'dopvar') % both are not dopvar type variables
%             Pcat = [a b];
%         elseif ~isa(a,'opvar')&&b.dim(1,1)==0 % a() is from R to L2, Note: L2 to L2 needs to be an opvar
%             if size(a,1)~=b.dim(2,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end 
% %             Pcat = b;
%             Pcat.Q2 = [a b.Q2];
%             Pcat.P = [zeros(b.dim(1,1),size(a,2)) b.P];
%         elseif ~isa(a,'opvar')&&b.dim(2,1)==0 % a() is from R to R, Note: L2 to R needs to be an opvar
%             if size(a,1)~=b.dim(1,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
% %             Pcat = b;
%             Pcat.P = [a b.P];
%             Pcat.Q2 = [zeros(b.dim(2,1),size(a,2)) b.Q2];
%         else
%         if any(b.dim(:,1)~=a.dim(:,1))
%             error('Cannot concatentate horizontally. A and B have different output dimensions');
%         end
% %         Pcat = b;
%         fset = {'P', 'Q1', 'Q2'};
% 
%         for i=fset
%             Pcat.(i{:}) = [a.(i{:}) b.(i{:})];
%         end
%         fset = {'R0','R1','R2'};
%         for i=fset
%             Pcat.R.(i{:}) = [a.R.(i{:}) b.R.(i{:})];
%         end
%         end
%     elseif ~isa(b,'dopvar')
%         if ~isa(b,'opvar')&&a.dim(1,1)==0 % b() is from R to L2, Note: L2 to L2 needs to be an opvar
%             if size(b,1)~=a.dim(2,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
% %             Pcat = a;
%             Pcat.Q2 = [a.Q2 b];
%             Pcat.P = [a.P zeros(a.dim(1,1),size(b,2))];
%         elseif ~isa(b,'opvar')&&a.dim(2,1)==0 % b() is from R to R, Note: L2 to R needs to be an opvar
%             if size(b,1)~=a.dim(1,1) 
%                 error('Cannot concatentate horizontally. A and B have different output dimensions');
%             end
% %             Pcat = a;
%             Pcat.P = [a.P b];
%             Pcat.Q2 = [a.Q2 zeros(a.dim(2,1),size(b,2))];
%         else
%         if any(b.dim(:,1)~=a.dim(:,1))
%             error('Cannot concatentate horizontally. A and B have different output dimensions');
%         end
% %         Pcat = a;
%         fset = {'P', 'Q1', 'Q2'};
% 
%         for i=fset
%             Pcat.(i{:}) = [a.(i{:}) b.(i{:})];
%         end
%         fset = {'R0','R1','R2'};
%         for i=fset
%             Pcat.R.(i{:}) = [a.R.(i{:}) b.R.(i{:})];
%         end
%         end
%     else
%         if any(b.dim(:,1)~=a.dim(:,1))
%             error('Cannot concatentate horizontally. A and B have different output dimensions');
%         end
% %         Pcat = a;
%         fset = {'P', 'Q1', 'Q2'};
% 
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
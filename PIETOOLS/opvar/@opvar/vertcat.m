function [Pcat] = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = vertcat(varargin) takes n inputs and concatentates them vertically,
% provided they satisfy the following criteria:
% 1) At least one input must be of type 'opvar', remaining inputs can be of
%    type 'double', 'polynomial', or 'opvar'.
% 2) The input dimensions varargin{j}.dim(:,2) of all opvar objects must
%       match;
% 3) The spatial variables varargin{j}.var1 and varargin{j}.var2, as well 
%       as the domain varargin{j}.I of all opvar objects must match;
% 4) Concatenation should make sense within the context of the opvar class,
%       that is, opvars always mapt to from RxL2, never to L2xR. We cannot
%       concatenate e.g. [A;B] for A:R-->L2 and B:R-->R.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - vertcat
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
% SS - 10/09/24: Revert to allow some matrix-opvar concatenations.
% DB - 11/13/24: "ambiguous" concatenations allowed, warning displayed.

% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% Extract the operators
a = varargin{1};    b = varargin{2};

% Currently support some matrix-opvar concatenation
if ~isa(b,'opvar')
    % Only supported if a in fact corresponds to a matrix as well
    if all(a.dim(2,:)==0)
        opvar tmp; tmp.I = a.I; tmp.var1 = a.var1; tmp.var2 = a.var2;
        tmp.P = b;
        b = tmp;
    else
        error("Ambiguous concatenation [a;b] of opvars. Please explicitly define 'a' and 'b' as opvars to resolve.");
    end
elseif ~isa(a,'opvar') && (isa(a,'double')||isa(a,'polynomial')||isa(a,'dpvar'))
    opvar tmp; tmp.I = b.I; tmp.var1 = b.var1; tmp.var2 = b.var2;
    if b.dim(2,2)~=0 % b from L2, so [a;b] is possible only if a is integral Q1 (or multiplier R0, but we ignore that)
        tmp.Q1 = a;
    elseif b.dim(1,2)~=0 && (isa(a,'double')||isempty(a.degmat)) % b maps from R, [a;b] is possible only if a is multiplier P (or Q2, but we ignore that)
        tmp.P = a;
    else
        error("Ambiguous concatenation [a;b] of opvars. Explicitly define 'a' and 'b' as opvars to resolve.");
    end
    tmp.dim = tmp.dim; % rectify dimensions
    a = tmp;
    clear tmp;
elseif ~isa(a,'opvar')
    error("Concatenation of 'opvar' objects is only supported with objects of type 'opvar', 'polynomial', or 'double'.")
end


% Check that domain and variables match
if any(a.I~=b.I)|| ~strcmp(a.var1.varname,b.var1.varname) || ~strcmp(a.var2.varname,b.var2.varname)
    error('Operators being concatenated have different intervals or different independent variables');
end
% Check that the output dimensions match
a.dim = a.dim;  b.dim = b.dim;  % correction to make components have consistent dimensions 8/27-ss
if any(a.dim(:,2)~=b.dim(:,2))
    error('Cannot concatenate vertically: Input dimensions of opvar objects do not match')
end

% Avoid "ambiguous" concatenations
if (a.dim(2,1)~=0 && b.dim(1,1)~=0)
    % a has rows mapping to L2, but b has rows mapping to R.
    % Concatenation would place those rows of a below those of b in the
    % opvar, which we currently prohibit...
    warning('Proposed opvar concatenation is ambiguous.')
end

% Finally, let's actually concatenate
opvar Pcat;
Pcat.I = a.I;   Pcat.var1 = a.var1;     Pcat.var2 = a.var2;
fset = {'P', 'Q1', 'Q2'};
for i=fset
    Pcat.(i{:}) = [a.(i{:}); b.(i{:})];
end
fset = {'R0','R1','R2'};
for i=fset
    Pcat.R.(i{:}) = [a.R.(i{:}); b.R.(i{:})];
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
% provided they satisfy the following criteria:
% 1) Atleast one input is an opvar variable.
% 2) If all the inputs are not opvar, then the operator maps from R to
% RxL2 or L2 to L2. 
% 3) Currently, it supports RxL2 to RxL2 concatenation only if ALL the inputs are
% opvar.
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
%         Pcat = vertcat(a,b);
%     elseif isa(b,'dpvar')
%         a = opvar2dopvar(a);
%         Pcat = vertcat(a,b);
%     elseif ~isa(a,'opvar') 
%         if ~isa(b,'opvar')
%             if size(a,2)~=size(b,2)
%                 error('Cannot concatentate vertically. A and B have different input dimensions');
%             end
%             Pcat = [a;b];
%         else
%             bdim = b.dim;
%             if size(a,2)~=sum(bdim(:,2))
%                 error("Cannot concatentate vertically. A and B have different input dimensions");
%             end
% %             Pcat = b;
%             if bdim(2,2) ==0
%                 % b:R-->RxL2, assume a:R-->R
%                 if isa(a,'polynomial') || isa(b,'dpvar')
%                     try a=double(a);
%                     catch
%                         error('Convert arguments to opvar for this concatenation')
%                     end
%                 end
%                 Pcat.P = [a; b.P]; % a() is from R to R
%                 Pcat.Q1 = [zeros(size(a,1),b.dim(2,2)); b.Q1];
%                 Pcat.Q2 = b.Q2;
%             elseif bdim(1,2)==0 
%                 % b:L2-->RxL2, assume a:L2-->R
%                 % NOTE: a:L2-->L2 must be specified as opvar
%                 %Pcat.Q2 = [zeros(size(a,1),bdim(1,2)); b.Q2];
%                 Pcat.Q1 = [a; b.Q1];
%                 Pcat.R.R0 = b.R.R0;
%                 Pcat.R.R1 = b.R.R1;
%                 Pcat.R.R2 = b.R.R2;
% %                 Pcat.Q1 = b.Q1;
% %                 Pcat.R.R0 = [a; b.R.R0];
% %                 Pcat.R.R1 = [zeros(size(a)); b.R.R1];
% %                 Pcat.R.R2 = [zeros(size(a)); b.R.R2];
%             else %find if such a operation is valid is any useful scenario and implement it
%                 error('Cannot concatenate vertically. This feature is not yet supported.');
%             end
%         end
%     elseif ~isa(b,'opvar') 
%         adim = a.dim;
%         if size(b,2)~=sum(adim(:,2))
%             error("Cannot concatentate vertically. A and B have different input dimensions");
%         end
% %         Pcat = a;
%         if adim(2,2) ==0
%             % a:R-->RxL2, assume b:R-->R
%             if isa(b,'polynomial') || isa(b,'dpvar')
%                 try b=double(b);
%                 catch
%                     error('Convert arguments to opvar for this concatenation')
%                 end
%             end
%             Pcat.P = [a.P; b]; % b() is from R to R
%             Pcat.Q1 = [a.Q1; zeros(size(b,1),a.dim(2,2))];
%             Pcat.Q2 = a.Q2;
%         elseif adim(1,2)==0 
%             % a:L2-->RxL2, assume b:L2-->R
%             % NOTE: b:L2-->L2 must be specified as opvar
%             %Pcat.Q2 = [a.Q2; zeros(size(b,1),adim(1,2))];
%             Pcat.Q1 = [a.Q1; b];
%             Pcat.R.R0 = a.R.R0;
%             Pcat.R.R1 = a.R.R1;
%             Pcat.R.R2 = a.R.R2;
% %             Pcat.Q1 = a.Q1;
% %             Pcat.R.R0 = [a.R.R0; b];
% %             Pcat.R.R1 = [a.R.R1; zeros(size(b))];
% %             Pcat.R.R2 = [a.R.R2; zeros(size(b))];
%         else %find if such a operation is valid is any useful scenario and implement it
%             error('Cannot concatenate vertically. This feature is not yet supported.');
%         end
%     else
%         if any(b.dim(:,2)~=a.dim(:,2))
%             error("Cannot concatentate vertically. A and B have different input dimensions");
%         end
% %         Pcat = a;
%         fset = {'P', 'Q1', 'Q2'};
%         for i=fset
%             Pcat.(i{:}) = [a.(i{:}); b.(i{:})];
%         end
%         fset = {'R0','R1','R2'};
%         for i=fset
%             Pcat.R.(i{:}) = [a.R.(i{:}); b.R.(i{:})];
%         end
%     end
%     if nargin>2 % Continue concatenation if inputs are more than 2
%         Pcat = vertcat(Pcat, varargin{3:end});
%     end
% end
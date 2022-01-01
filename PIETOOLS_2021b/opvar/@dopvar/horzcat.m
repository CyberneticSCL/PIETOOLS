function [Pcat] = horzcat(varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = horzcat(varargin) takes n-inputs and concatentates them horizontally,
% provided they satisfy the following criterias.
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


if nargin==1
    Pcat = varargin{1};
else
    a = varargin{1};
    b = varargin{2};
    
    if isa(a,'dopvar') % correction to make components have consistent dimensions 8/27-ss
        a.dim = a.dim;
    end
    if isa(b,'dopvar') % correction to make components have consistent dimensions 8/27-ss
        b.dim = b.dim;
    end
    
    
    
    if ~isa(a,'dopvar')
        if ~isa(b,'dopvar') % both are not dopvar type variables
            Pcat = [a b];
        elseif ~isa(a,'opvar')&&b.dim(1,1)==0 % a() is from R to L2, Note: L2 to L2 needs to be an opvar
            if size(a,1)~=b.dim(2,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end 
            Pcat = b;
            Pcat.Q2 = [a b.Q2];
            Pcat.P = [zeros(b.dim(1,1),size(a,2)) b.P];
        elseif ~isa(a,'opvar')&&b.dim(2,1)==0 % a() is from R to R, Note: L2 to R needs to be an opvar
            if size(a,1)~=b.dim(1,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end
            Pcat = b;
            Pcat.P = [a b.P];
            Pcat.Q2 = [zeros(b.dim(2,1),size(a,2)) b.Q2];
        else
        if any(b.dim(:,1)~=a.dim(:,1))
            error('Cannot concatentate horizontally. A and B have different output dimensions');
        end
        Pcat = b;
        fset = {'P', 'Q1', 'Q2'};

        for i=fset
            Pcat.(i{:}) = [a.(i{:}) b.(i{:})];
        end
        fset = {'R0','R1','R2'};
        for i=fset
            Pcat.R.(i{:}) = [a.R.(i{:}) b.R.(i{:})];
        end
        end
    elseif ~isa(b,'dopvar')
        if ~isa(b,'opvar')&&a.dim(1,1)==0 % b() is from R to L2, Note: L2 to L2 needs to be an opvar
            if size(b,1)~=a.dim(2,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end
            Pcat = a;
            Pcat.Q2 = [a.Q2 b];
            Pcat.P = [a.P zeros(a.dim(1,1),size(b,2))];
        elseif ~isa(b,'opvar')&&a.dim(2,1)==0 % b() is from R to R, Note: L2 to R needs to be an opvar
            if size(b,1)~=a.dim(1,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end
            Pcat = a;
            Pcat.P = [a.P b];
            Pcat.Q2 = [a.Q2 zeros(a.dim(2,1),size(b,2))];
        else
        if any(b.dim(:,1)~=a.dim(:,1))
            error('Cannot concatentate horizontally. A and B have different output dimensions');
        end
        Pcat = a;
        fset = {'P', 'Q1', 'Q2'};

        for i=fset
            Pcat.(i{:}) = [a.(i{:}) b.(i{:})];
        end
        fset = {'R0','R1','R2'};
        for i=fset
            Pcat.R.(i{:}) = [a.R.(i{:}) b.R.(i{:})];
        end
        end
    else
        if any(b.dim(:,1)~=a.dim(:,1))
            error('Cannot concatentate horizontally. A and B have different output dimensions');
        end
        Pcat = a;
        fset = {'P', 'Q1', 'Q2'};

        for i=fset
            Pcat.(i{:}) = [a.(i{:}) b.(i{:})];
        end
        fset = {'R0','R1','R2'};
        for i=fset
            Pcat.R.(i{:}) = [a.R.(i{:}) b.R.(i{:})];
        end
    end
    if nargin>2 % check if there are more than 2 objects that need to stacked
        Pcat = horzcat(Pcat, varargin{3:end});
    end
end
end
function [Pcat] = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = vertcat(varargin) takes n-inputs and concatentates them vertically,
% provided they satisfy the following criterias.
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



if nargin==1
    Pcat = varargin{1};
else
    a = varargin{1};
    b = varargin{2};
    
    if isa(a,'opvar') % correction to make components have consistent dimensions 8/27-ss
        a.dim = a.dim;
    end
    if isa(b,'opvar') % correction to make components have consistent dimensions 8/27-ss
        b.dim = b.dim;
    end
    
    if ~isa(a,'opvar') 
        if ~isa(b,'opvar')
            if size(a,2)~=size(b,2)
                error('Cannot concatentate vertically. A and B have different input dimensions');
            end
            Pcat = [a;b];
        else
            bdim = b.dim;
            if size(a,2)~=sum(bdim(:,2))
                error("Cannot concatentate vertically. A and B have different input dimensions");
            end
            Pcat = b;
            if bdim(2,2) ==0
                Pcat.P = [a; b.P]; % a() is from R to R
                Pcat.Q1 =[zeros(size(a,1),b.dim(2,2)); b.Q1];
            elseif bdim(1,2)==0 % a() is from L2 to L2, Note: a() cannot be a matrix and map L2 to R 
                Pcat.Q2 = [zeros(size(a,1),b.dim(1,2)); b.Q2];
                Pcat.R.R0 = [a; b.R.R0];
                Pcat.R.R1 = [zeros(size(a)); b.R.R1];
                Pcat.R.R2 = [zeros(size(a)); b.R.R2];
            else %find if such a operation is valid is any useful scenario and implement it
                error('Cannot concatenate vertically. This feature is not yet supported.');
            end
        end
    elseif ~isa(b,'opvar') 
        adim = a.dim;
        if size(b,2)~=sum(adim(:,2))
            error("Cannot concatentate vertically. A and B have different input dimensions");
        end
        Pcat = a;
        if adim(2,2) ==0
            Pcat.P = [a.P; b]; % b() is from R to R
            Pcat.Q1 = [a.Q1; zeros(size(b,1),a.dim(2,2))];
        elseif adim(1,2)==0 % b() is from L2 to L2, Note: b() cannot be L2 to R and not be opvar
            Pcat.Q2 = [a.Q2; zeros(size(b,1),a.dim(1,2))];
            Pcat.R.R0 = [a.R.R0; b];
            Pcat.R.R1 = [a.R.R1; zeros(size(b))];
            Pcat.R.R2 = [a.R.R2; zeros(size(b))];
        else %find if such a operation is valid is any useful scenario and implement it
            error('Cannot concatenate vertically. This feature is not yet supported.');
        end
    else
        if any(b.dim(:,2)~=a.dim(:,2))
            error("Cannot concatentate vertically. A and B have different input dimensions");
        end
        Pcat = a;
        fset = {'P', 'Q1', 'Q2'};
        for i=fset
            Pcat.(i{:}) = [a.(i{:}); b.(i{:})];
        end
        fset = {'R0','R1','R2'};
        for i=fset
            Pcat.R.(i{:}) = [a.R.(i{:}); b.R.(i{:})];
        end
    end
    if nargin>2 % Continue concatenation if inputs are more than 2
        Pcat = vertcat(Pcat, varargin{3:end});
    end
end
end
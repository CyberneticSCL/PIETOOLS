function varargout = size(T,dim1,dim2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varargout = size(T,dim1,dim2) determines the dimensions of an opvar2d object 
%
% Version 1.0
% Date: 07/12/21
% 
% INPUT
% T:    dopvar2d class objects for which to establish the dimensions
% dim:  (optional) inputs indicating the dimension along which the size
%       is desired
%
% OUTPUT
% The size of the object along the specified dimensions
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - size
%
% Copyright (C)2021  M. Peet, S. Shivakumar
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
% Initial coding DJ - 07_12_2021

% For only one input, return the full number of rows and columns of the
% dopvar2d object
if nargin ==1
    Tdim = T.dim;
    out = [sum(Tdim(:,1)),sum(Tdim(:,2))];
    if nargout==0 || nargout==1
        varargout{1} = out;
    elseif nargout==2
        varargout{1} = out(1);
        varargout{2} = out(2);
    end  
% For two inputs, return either the full number of rows or columns of the
% opvar2d object, or the full dimensions of the object
elseif nargin == 2
    if strcmp(dim1,':') || strcmp(dim1,'all')
        out = T.dim;
        if nargout==0 || nargout==1
            varargout{1} = out;
        elseif nargout==2
            varargout{1} = out(:,1);
            varargout{2} = out(:,2);
        elseif nargout==4
            varargout{1} = out(1,:);
            varargout{2} = out(2,:);
            varargout{3} = out(3,:);
            varargout{4} = out(4,:);
        end  
    elseif ~( isa(dim1,'double') && isreal(dim1) &&  all(size(dim1)==[1 1]) &&...
            all(ceil(dim1)==floor(dim1)) )
        error('Dimension must be a positive integer scalar');
    elseif all(size(dim1)==[1 1])
        if dim1==1 || dim1==2
            Tdim = T.dim;
            varargout{1} = sum(Tdim(:,dim1));
        elseif dim1>0
            varargout{1} = 1;
        else
            error('Dimension must be a positive integer scalar');
        end
    end
% For three inputs, return a slice of the full dimensions of the dopvar2d
% object
elseif nargin == 3
    if strcmp(dim1,':') || strcmp(dim1,'all')
        if strcmp(dim2,':') || strcmp(dim2,'all')
            varargout = size(T,':');
        elseif ~( isa(dim2,'double') && isreal(dim2) &&  all(size(dim2)==[1 1]) &&...
                all(ceil(dim2)==floor(dim2)) )
            error('Second dimension must be a positive integer scalar');
        elseif all(size(dim2)==[1 1])
            if dim2==1 || dim2==2
                Tdim = T.dim;
                varargout{1} = Tdim(:,dim2);
            elseif dim2>0
                varargout{1} = 1;
            else
                error('Second dimension must be a positive integer scalar');
            end
        end
    elseif ~( isa(dim1,'double') && isreal(dim1) &&  all(size(dim1)==[1 1]) &&...
            all(ceil(dim1)==floor(dim1)) )
        error('First dimension must be a positive integer scalar');
    elseif strcmp(dim2,':') || strcmp(dim2,'all')
        if dim1==1 || dim1==2 || dim1==3 || dim1==4
            Tdim = T.dim;
            varargout{1} = Tdim(dim1,:);
        else
            error('First dimension must be one of [1,2,3,4]');
        end
    elseif ~( isa(dim2,'double') && isreal(dim2) &&  all(size(dim2)==[1 1]) &&...
                all(ceil(dim2)==floor(dim2)) )
        error('Dimension must be a positive integer scalar');
    elseif dim1==1 || dim1==2 || dim1==3 || dim1==4
        if dim2==1 || dim2==2
            Tdim = T.dim;
            varargout{1} = Tdim(dim1,dim2);
        elseif dim2>0
            varargout{1} = 1;
        else
            error('Second dimension must be a positive integer scalar');
        end
    else
        error('First dimension must be one of [1,2,3,4]');
    end
end

end  
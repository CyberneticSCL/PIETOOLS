function [params,out] = convert(obj,convertTo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that converts obj from one type to another (typically
% PIE)
% Input: 
% obj - sys class object
% convertTo - string, 'pie' or 'ddf'
% Output:
% obj - sys class object with pie or ddf parameters
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
arguments
    obj;
    convertTo {mustBeMember(convertTo,{'pie','ddf'})} = '';
end

if strcmp(convertTo, 'pie')
    if strcmp(obj.type,'pde')
        tmp = convert_PIETOOLS_PDE(obj.params);
    elseif strcmp(obj.type,'dde')
        tmp = convert_PIETOOLS_DDE(obj.params,'pie');
    elseif strcmp(obj.type,'nds')
        tmp = convert_PIETOOLS_NDS(obj.params,'pie');
    elseif strcmp(obj.type,'ddf')
        tmp = convert_PIETOOLS_DDF(obj.params,'pie');
    elseif strcmp(obj.type,'pie')
        out = obj;
        return;
    end
    out = sys();
    out.type = 'pie';
    out.params = tmp;
    params = tmp;
elseif strcmp(convertTo,'ddf') && (strcmp(obj.type,'nds')||strcmp(obj.type,'dde'))
    out = obj;
    if strcmp(obj.type,'nds')
        tmp = convert_PIETOOLS_NDS2DDF(obj.params);
    elseif strcmp(obj.type,'dde')
        tmp = minimize_PIETOOLS_DDE2DDF(obj.params);
    end
    out.type = 'ddf';
    out.params = tmp;
    params = tmp;
end

fprintf('Conversion to %s was successful\n', convertTo);
end
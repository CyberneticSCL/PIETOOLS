function out=getprop(T,prop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out=get(T,prop) retrieves specified components of the dopvar2d object T
% 
% INPUTS:
% T : A dopvar2d variable
% prop: structure with property to be accessed
% 
% OUTPUTS:
% out: property value
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if length(prop)>5
    error('Invalid property name: No more than 5 layers of properties can be extracted');
end

% The first property must refer to one of the fields of the structure (R00
% up to R22)
switch prop(1).type
    case '.'
        prop1 = prop(1).subs;
        if strcmp(prop1,'fieldnames') || strcmp(prop1,'fnames')
            out = {'I','dim','var1','var2','R00', 'R0x', 'R0y', 'R02', 'Rx0', 'Rxx', 'Rxy', 'Rx2', 'Ry0', 'Ryx', 'Ryy', 'Ry2', 'R20', 'R2x', 'R2y', 'R22'};
        elseif strcmp(prop1,'parameters') || strcmp(prop1,'params')
            out = {'R00', 'R0x', 'R0y', 'R02', 'Rx0', 'Rxx', 'Rxy', 'Rx2', 'Ry0', 'Ryx', 'Ryy', 'Ry2', 'R20', 'R2x', 'R2y', 'R22'};
        elseif strcmp(prop1,'oparameters') || strcmp(prop1,'oparams')
            out = {'R00', 'R0x', 'R0y', 'R02', 'Rx0', 'Rxy', 'Ry0', 'Ryx', 'R20'};
        elseif strcmp(prop1,'xparameters') || strcmp(prop1,'xparams')
            out = {'Rxx', 'Rx2', 'R2x'};
        elseif strcmp(prop1,'yparameters') || strcmp(prop1,'yparams')
            out = {'Ryy', 'Ry2', 'R2y'};
        elseif strcmp(prop1,'xyparameters') || strcmp(prop1,'xyparams')
            out = {'R22'};
        elseif strcmp(prop1,'dvarname') || strcmp(prop1,'dvarnames') || strcmp(prop1,'dvars')
            out = {};
            fset0 = {'R00', 'R0x', 'R0y', 'R02', 'Rx0', 'Rxy', 'Ry0', 'Ryx', 'R20'};
            for f=fset0
                tmp_R = T.(f{:});
                if ~isempty(tmp_R)
                    out = union(out,tmp_R.dvarname);
                end
            end
            fsetx = {'Rxx', 'Rx2', 'R2x'};
            for f=fsetx
                tmp_R = T.(f{:});
                for i=1:3
                    if isa(tmp_R{i},'dpvar')
                        out = union(out,tmp_R{i}.dvarname);
                    end
                end
            end
            fsety = {'Ryy', 'Ry2', 'R2y'};
            for f=fsety
                tmp_R = T.(f{:});
                for j=1:3
                    if isa(tmp_R{j},'dpvar')
                        out = union(out,tmp_R{j}.dvarname);
                    end
                end
            end
            fset2 = {'R22'};
            for f=fset2
                tmp_R = T.(f{:});
                for i=1:3
                    for j=1:3
                        if isa(tmp_R{i,j},'dpvar')
                            out = union(out,tmp_R{i,j}.dvarname);
                        end
                    end
                end
            end            
        else
            out = T.(prop1);
        end
    otherwise
        error('Invalid property name: First property must be one of the field names of the object');
end

% Remaining components may take cell elements, array elements, or fields of
% the desired object
if length(prop) >= 2
    for i=2:length(prop)
        switch prop(i).type
            case '.'
                propi = prop(i).subs;
                out = out.(propi);
            case '()'
                if length(prop(i).subs)==1  % use linear indexing
                    if strcmp(prop(i).subs{1},':')
                        indr = 1:size(out,1);   indc = 1:size(out,2);
                    else
                        [nr,nc] = size(out);
                        [indr,indc] = ind2sub([nr,nc],prop(i).subs{1});
                    end
                else
                    indr = prop(i).subs{1};
                    indc = prop(i).subs{2};
                end
                out = out(indr,indc);
            case '{}'
                if length(prop(i).subs)==1  % use linear indexing
                    if strcmp(prop(i).subs{1},':')
                        indr = 1:size(out,1);   indc = 1:size(out,2);
                    else
                        [nr,nc] = size(out);
                        [indr,indc] = ind2sub([nr,nc],prop(i).subs{1});
                    end
                else
                    indr = prop(i).subs{1};
                    indc = prop(i).subs{2};
                end
                out = out{indr,indc};
            otherwise
                error('Invalid property name: only types ''.'', ''()'', and ''{}'' are supported');
        end
    end
end
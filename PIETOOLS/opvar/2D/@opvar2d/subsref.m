function Psop=subsref(Pbop,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psop = subsref(Pbop,ref) extracts components or certain rows or columns
% of an opvar2d object
%
% Version 1.0
% Date: 07/05/21
% 
% INPUT
% Pbop: opvar class object to slice
% ref:  specification of the components to extract
% 
% OUTPUT
% Psop: slice of the opvar object
%  
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - subsref
%
% Copyright (C)2022 M. Peet, S. Shivakumar, D. Jagt
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07_05_2021
%   ^ Based heavily on "@opvar"-subsref code by SS ^
% DJ, 05/24/22: Allow domain to be extracted with dom;
% DJ, 08/07/22: Allow logical indexing;
%               Also allow parameter indexing using '({})';
% DJ, 10/20/2024: Add support for extracting variables as Pop.vars;

switch ref(1).type
    case '.'
        if strcmp(ref(1).subs,'vars')
            % Allow spatial variables to be extracted as Pbop.vars.
            var1 = Pbop.var1;
            var2 = Pbop.var2;
            vars = [var1,var2];
            if numel(ref)==1
                % Just Pbop.vars.
                Psop = vars;
            else
                % e.g. Pbop.vars(1).
                Psop = builtin('subsref',vars,ref(2:end));
            end
            return
        end
        if numel(ref)==1 && strcmp(ref.subs,'dom')
            % Allow spatial domain to be extracted as Pbop.dom;
            ref.subs = 'I';
        end
        Psop = getprop(Pbop,ref);
        
    case '()'
        dim = Pbop.dim;
        nr = sum(dim(:,1));     nc = sum(dim(:,2));
        if length(ref(1).subs)==1  
            % Convert linear indexing to row/column indexing.
            % Does linear indexing here actually make sense???
            if strcmp(ref(1).subs{1},':')
                indr = 1:dim(1);
                indc = 1:dim(2);
            elseif isa(ref(1).subs{1},'cell')
                error('Parameter indexing using linear indices is not supported.')
            else
                [indr,indc] = ind2sub([nr,nc],ref(1).subs{1});
            end
        elseif length(ref(1).subs)==2
            % Convert row/column indices to a standard format.
            if strcmp(ref(1).subs{1},':')
                indr = 1:nr;
            elseif isa(ref(1).subs{1},'cell')
                % Parameter indexing: extract parameters
                % R0., Rx., Ry., and/or R2..
                indr_param = cell2mat(ref(1).subs{1});
                if islogical(indr_param)
                    if length(indr_param)~=4
                        error('Logical parameter indexing requires 4 row indices.')
                    end
                else
                    if max(indr_param)>4
                        error('The proposed parameter row-indices exceed the number of output spaces (4).');
                    end
                end
                % Convert the parameter indices to standard row indices.
                nnr_op = cumsum([0;dim(:,1)]);
                nr_indcs = mat2cell([nnr_op(1:end-1),nnr_op(2:end)],ones(4,1));
                nr_indcs = cellfun(@(x) (x(1)+1:x(2))', nr_indcs, 'UniformOutput',false);
                nr_indcs = nr_indcs(indr_param);
                indr = cell2mat(nr_indcs)';
            elseif islogical(ref(1).subs{1})
                % Convert logical indices to standard indices.
                if length(ref(1).subs{1})~=nr
                    error('The number of logical row indices should match the total number of rows in the operator.')
                end
                indr = 1:nr;
                indr = indr(ref(1).subs{1});                    
            else
                % Extract standard indices.
                indr = ref(1).subs{1};
                if max(indr)>nr
                    error('The proposed row indices exceed the row dimension of the input object.');
                end
            end
            
            if strcmp(ref(1).subs{2},':')
                indc = 1:nc;
            elseif isa(ref(1).subs{2},'cell')
                % Parameter indexing: extract parameters
                % R.0., R.x, R.y, and/or R.2.
                indc_param = cell2mat(ref(1).subs{2});
                if islogical(indc_param)
                    if length(indc_param)~=4
                        error('Logical parameter indexing requires 4 column indices.')
                    end
                else
                    if max(indc_param)>4
                        error('The proposed parameter column-indices exceed the number of input spaces (4).');
                    end
                end
                % Convert the parameter indices to standard column indices.
                nnc_op = cumsum([0;dim(:,2)]);
                nc_indcs = mat2cell([nnc_op(1:end-1),nnc_op(2:end)],ones(4,1));
                nc_indcs = cellfun(@(x) (x(1)+1:x(2))', nc_indcs, 'UniformOutput',false);
                nc_indcs = nc_indcs(indc_param);
                indc = cell2mat(nc_indcs)';
            elseif islogical(ref(1).subs{2})
                % Convert logical indices to standard indices.
                if length(ref(1).subs{2})~=nc
                    error('The number of logical column indices should match the total number of columns in the operator.')
                end
                indc = 1:nc;
                indc = indc(ref(1).subs{2}); 
            else
                % Extract standard indices.
                indc = ref(1).subs{2};
            end
        else
            error('Opvar2d objects only take row and column indices.');
        end
        % Outsource actual extracting of the desired elements to op_slice.
        Psop = op_slice(Pbop,indr,indc);
        
    case '{}'
%         % Extract the parameter indexed by ref.subs.
%         Rparams = {'R00', 'R0x', 'R0y', 'R02';
%                    'Rx0', 'Rxx', 'Rxy', 'Rx2';
%                    'Ry0', 'Ryx', 'Ryy', 'Ry2';
%                    'R20', 'R2x', 'R2y', 'R22'};
%         if length(ref(1).subs)==1
%             if islogical(ref(1).subs{1})
%                 error('Linear indexing using logical values is not supported.')
%             else
%                 if any(ref(1).subs{1}>16)
%                     error('For cell indexing, linear indices must range between 1 and 16.')
%                 end
%                 indp = ref(1).subs{1};
%             end
%         elseif length(ref(1).subs)==2
%             % Convert to linear indices.
%             if islogical(ref(1).subs{1})
%                 if length(ref(1).subs{1})~=4
%                     error('For cell indexing, the number of logical row indices must be equal to 4.')
%                 end
%                 indr = 1:4;
%                 indr = indr(ref(1).subs{1});
%             else
%                 if any(ref(1).subs{1}>4)
%                     error('For cell indexing, row and column indices must range between 1 and 4.')
%                 end
%                 indr = ref(1).subs{1};
%             end
%             if islogical(ref(1).subs{2})
%                 if length(ref(1).subs{2})~=4
%                     error('For cell indexing, the number of logical column indices must be equal to 4.')
%                 end
%                 indc = 1:4;
%                 indc = indc(ref(1).subs{2});
%             else
%                 if any(ref(1).subs{2}>4)
%                     error('For cell indexing, row and column indices must range between 1 and 4.')
%                 end
%                 indc = ref(1).subs{2};
%             end
%             if length(indr)~=length(indc)
%                 error('The number of row and column indices should match.')
%             end
%             indp = indr + (indc-1)*4;
%         else
%             error('Only row and column indices are supported for opvar2d objects.')
%         end
        
        error('Indexing using ''{}'' is currently not supported');
end
end
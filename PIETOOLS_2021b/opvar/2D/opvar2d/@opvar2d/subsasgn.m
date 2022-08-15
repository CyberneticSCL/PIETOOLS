function b = subsasgn(a,L,RHS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b = subsasgn(a,L,RHS) assigns values RHS to field/slice L of an opvar2d 
% object a
%
% Version 1.0
% Date: 07/06/21
% 
% INPUT
% a:    opvar2d class object
% L:    a struct specifying the component to be adjusted/assigned
% RHS:  a value to be assigned to the component
%
% OUTPUT
% b:    opvar2d object with desired component value
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or D. Jagt at djagt@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07_06_2021
% DJ, 08/07/2022: Allow '()'-type subsref;
%                 Also allow parameter indexing using '({})';

%a = opvar2d(a);
%sza = size(a);
switch L(1).type
    case '.'
        if length(L) == 1
            temp = RHS;
        else
            % Peform all subsasgn but L(1)
            temp = subsref(a,L(1));
            switch L(2).type
                case '{}'
                    if length(L(2).subs)==1
                        [nr,nc] = size(temp);
                        [rindx,cindx] = ind2sub([nr,nc],L(2).subs{1});
                    elseif length(L(2).subs)>=2
                        rindx = L(2).subs{1};   cindx = L(2).subs{2};
                    end
                    if length(L)==2
                        temp{rindx,cindx} = RHS;
                    else
                        temp{rindx,cindx} = subsasgn(temp{rindx,cindx},L(3:end),RHS);
                    end
                otherwise
                    temp = subsasgn(temp,L(2:end),RHS);
            end
        end
        b = set(a,L(1).subs,temp);
        
    case '()'
        if ~isa(RHS,'opvar2d')
            error('Type ''()'' subsasgn requires the right-hand side to be an opvar2d object')
        else
            % Make sure we can work with the inputs.
            if ~isequal(a.I,RHS.I) || ~all(isequal(a.var1,RHS.var1)) || ~all(isequal(a.var2,RHS.var2))
                error('The spatial domain and variables of the right-hand side operator should match those of the old operator.')
            elseif length(L(1).subs)==1
                error('Type ''()'' subsasgn using linear indices is not supported.')
            elseif length(L(1).subs)>=3
                error('Opvar2d objects take at most two subscripts.')
            end
            
            % Extract the dimensions of the operators a and RHS.
            a_dim = a.dim;
            nr_a = sum(a_dim(:,1));     nc_a = sum(a_dim(:,2));
            nnr_op_a = cumsum([0;a_dim(:,1)]);
            nnc_op_a = cumsum([0;a_dim(:,2)]);   
            RHS_dim = RHS.dim;
            nr_RHS = sum(RHS_dim(:,1)); nc_RHS = sum(RHS_dim(:,2));
            
            % Initialize the desired row indices in a.
            indr = L(1).subs{1};
            indr_param = 1:4;  % Assume all parameters are adjusted for now.
            if strcmp(indr,':')
                indr = 1:nr_a;
            elseif isa(indr,'cell')
                % Parameter indexing: adjust parameters
                % R0., Rx., Ry., and/or R2..
                indr_param = cell2mat(L(1).subs{1});
                if islogical(indr_param)
                    % Convert logical indices to standard indices.
                    if length(indr_param)~=4
                        error('Logical parameter indexing requires 4 row indices.')
                    end
                    indr_param_new = 1:nr_a;
                    indr_param = indr_param_new(indr_param);
                else
                    if max(indr_param)>4
                        error('The proposed parameter row-indices exceed the number of output spaces (4).');
                    elseif length(unique(indr_param))~=length(indr_param)
                        error('The parameter row indices should be unique');
                    end
                    indr_param = sort(indr_param);  % the order of the parameters is fixed.
                end
                % Convert the parameter indices to standard row indices.
                nr_indcs = mat2cell([nnr_op_a(1:end-1),nnr_op_a(2:end)],ones(4,1));
                nr_indcs = cellfun(@(x) (x(1)+1:x(2))', nr_indcs, 'UniformOutput',false);
                nr_indcs = nr_indcs(indr_param);
                indr = cell2mat(nr_indcs)';
            elseif islogical(indr)
                % Convert logical indices to standard indices.
                if length(indr)~=nr_a
                    error('The number of logical row indices should match the row dimension of the input object.');
                end
                indr = 1:nr_a;
                indr = indr(L(1).subs{1});
            else
                if max(indr)>nr_a
                    error('The proposed row indices exceed the row dimension of the input object.');
                elseif length(unique(indr))~=length(indr)
                    error('The row indices should be unique');
                end
                indr = indr(:)';
            end
            if length(indr)~=nr_RHS
                error('The number of row indices does not match the row dimension of the proposed right-hand side.')
            end
            
            % Initialize the desired column indices in RHS.
            indc = L(1).subs{2};
            if strcmp(indr,':')
                indc = 1:nc_a;
            elseif isa(indc,'cell')
                % Parameter indexing: adjust parameters
                % R.0, R.x, R.y, and/or R.2.
                indc_param = cell2mat(L(1).subs{2});
                if islogical(indc_param)
                    % Convert logical indices to standard indices.
                    if length(indc_param)~=4
                        error('Logical parameter indexing requires 4 row indices.')
                    end
                    indc_param_new = 1:nc_a;
                    indc_param = indc_param_new(indc_param);
                else
                    if max(indc_param)>4
                        error('The proposed parameter column-indices exceed the number of output spaces (4).');
                    elseif length(unique(indc_param))~=length(indc_param)
                        error('The parameter column-indices should be unique.');
                    end
                    indc_param = sort(indc_param);  % the order of the parameters is fixed.
                end
                % Convert the parameter indices to standard column indices.
                nc_indcs = mat2cell([nnc_op_a(1:end-1),nnc_op_a(2:end)],ones(4,1));
                nc_indcs = cellfun(@(x) (x(1)+1:x(2))', nc_indcs, 'UniformOutput',false);
                nc_indcs = nc_indcs(indc_param);
                indc = cell2mat(nc_indcs)';
            elseif islogical(indc)
                % Convert logical indices to standard indices.
                if length(indc)~=nc_a
                    error('The number of logical column indices should match the column dimension of the input object.');
                end
                indc = 1:nc_a;
                indc = indc(L(1).subs{2});
            else
                if max(indc)>nc_a
                    error('The proposed column indices exceed the column dimension of the input object.');
                elseif length(unique(indc))~=length(indc)
                    error('The column indices should be unique');
                end
                indc = indc(:)';
            end
            if length(indc)~=nc_RHS
                error('The number of column indices does not match the column dimension of the proposed right-hand side.')
            end
            
        end
        % Establish which row and column indices correspond to which
        % parameters in the operators a and RHS.     
        nnr_op_RHS = cumsum([0;RHS_dim(:,1)]);
        nnc_op_RHS = cumsum([0;RHS_dim(:,2)]);
        
        nnr_op_a = nnr_op_a([indr_param(:);indr_param(end)+1]);
        nnc_op_a = nnc_op_a([indc_param(:);indc_param(end)+1]);
        nnr_op_RHS = nnr_op_RHS([indr_param(:);indr_param(end)+1]);
        nnc_op_RHS = nnc_op_RHS([indc_param(:);indc_param(end)+1]);
        
        % Loop over each of the parameters in the operator, adjusting the
        % appropriate elements of each parameter to match those of the
        % parameters in the RHS operator.
        Rparams = {'R00', 'R0x', 'R0y', 'R02';
                   'Rx0', 'Rxx', 'Rxy', 'Rx2';
                   'Ry0', 'Ryx', 'Ryy', 'Ry2';
                   'R20', 'R2x', 'R2y', 'R22'};
        Rparams = Rparams(indr_param,indc_param);   % exclude pararmeters we know are not adjusted.
        for ll=1:numel(Rparams)
            % Establish which of the proposed rows and columns of the
            % operator RHS appear in the parameter indexed by ll.
            [r_param, c_param] = ind2sub(size(Rparams),ll);
            RHS_rindcs = find(indr>nnr_op_a(r_param) & indr<=nnr_op_a(r_param+1));
            RHS_cindcs = find(indc>nnc_op_a(c_param) & indc<=nnc_op_a(c_param+1));
            
            if isempty(RHS_rindcs) || isempty(RHS_cindcs)
                % If no adjustments need to be made to the parameter
                % indexed by ll, move on to the next paramter.
                continue
            elseif any(RHS_rindcs<=nnr_op_RHS(r_param) | RHS_rindcs>nnr_op_RHS(r_param+1))
                % Make sure that the rows in a correspond to the same type
                % of parameter in RHS as well.
                error('The output spaces the proposed right-hand side maps to should match those of the associated rows in the operator')
            elseif any(RHS_cindcs<=nnc_op_RHS(c_param) | RHS_cindcs>nnc_op_RHS(c_param+1))
                error('The output spaces the proposed right-hand side maps to should match those of the associated rows in the operator')
            end
            % Establish which of the proposed rows and columns of the
            % operator a must be adjusted.
            a_rindcs = indr(RHS_rindcs);
            a_cindcs = indc(RHS_cindcs);
            
            % Adjust the indices to correspond only to those of the
            % parameter ll.
            RHS_rindcs = RHS_rindcs - nnr_op_RHS(r_param);
            RHS_cindcs = RHS_cindcs - nnc_op_RHS(c_param);
            a_rindcs = a_rindcs - nnr_op_a(r_param);
            a_cindcs = a_cindcs - nnc_op_a(c_param);
            
            if isa(a.(Rparams{ll}),'cell')
                for kk=1:numel(a.(Rparams{ll}))
                    a.(Rparams{ll}){kk}(a_rindcs,a_cindcs) = RHS.(Rparams{ll}){kk}(RHS_rindcs,RHS_cindcs);
                end                       
            else
                a.(Rparams{ll})(a_rindcs,a_cindcs) = RHS.(Rparams{ll})(RHS_rindcs,RHS_cindcs);
            end
        end
        % Set the new operator.
        b = a;
        
    case '{}'
        error('{}- like subsassign is not supported for opvar2d objects.');
end

end

classdef (InferiorClasses={?polynomial,?dpvar})sopvar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - sopvar class
% Represents PI maps from L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
%
% This defines PI operators from one L2 space to another. 
%
%   Pop: L_2^p[Si,...Sn] --> L_2^q[Sj,...,Sm]
%
% Elements of this class are intended to be included in a container object.
% 
%   The operator has the form 
%         sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum) Z_d(S_2,S_3) C Z_d(S_3dum,S_1)
%
%   where S_1 is the variables in S_i (output) not in S_j (input)
%         S_2 is the variables in S_j (input) not in S_i (output)
%         S_3 is the variables common to S_j (input) and S_i (output)
%         S_3dum are dummy versions of the variables in S_3
%
% CLASS properties
% - P.C:   3^N x 1 Cell array of Quadpoly objects, where N is the number of variables which appear in both the input and output space;
%           Each quadpoly has the structure (I X Z_d(S_2,S_3)) C (I X Z_d(S_3dum,S_1))
%           NOTE: Wouldnt it be better to make this a NxNxN array so we don't have to figure out where everything is? Not as fast maybe, but the bottleneck will be quadvar operations anyway. 
%                Sachin: It is not (3^N x 1) cell, it is 3x3x...x3 where
%                        some of the 3's are replaced by 1's when that particular
%                        variable is not present in both domain and range.
% - P.dims:   2x1 vector specifying the rows and columns of the operator i.e. [p,q];
% - P.vars_in: nx1 cell array specifying the unique names of the spatial variables which appear in the input space
% - P.vars_out: mx1 cell array specifying the unique names of the spatial variables which appear in the output space
% - P.dom_in: nx2 array specifying the spatial domain of each variable in input space, so that dom_in(i,:) = [ai,bi]
% - P.dom_out: mx2 array specifying the spatial domain of each variable in output space, so that dom_in(i,:) = [ai,bi]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2026 PIETOOLS Team
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
% SS,MP, 01/15/2026: Initial coding
    properties
        vars_in = {};
        dom_in = zeros(0, 2);
        vars_out = {};
        dom_out = zeros(0, 2);
        dims = [1,1];
        params = {quadPoly(0,[],[],[1,1],{},{})};
        % each value in the cell is a quadPoly
        % for maps to and from \R spaces, add an extra dimension?
    end

    methods
        function P = sopvar(varsin,varsout,dims,dom_in,dom_out,params)
            P.vars_in = varsin(:).';
            P.vars_out = varsout(:).';
            P.dims = dims;
            P.dom_in=dom_in;
            P.dom_out=dom_out;
            if nargin==5                                                    % create a 0-valued operator with specified domains

                vars = union(P.vars_out, P.vars_in);                        % list of unique spatial variables
                n= numel(vars);                                             % total number of unique spatial variables
                celldim = repmat(3,1,n);
                [~,idx] = ismember(vars, intersect(P.vars_in,P.vars_out));  % position of common input-output variables in vars
                celldim(~idx) = 1;
                if n == 0
                    P.params = repmat({zeros(dims(1),dims(2))},1,1);
                else
                    P.params = repmat({zeros(dims(1),dims(2))},celldim);
                end
            else
                P.params = params;
            end
        end
    end
    
    methods(Static)
        out = randOpvar(varsin,varsout,dim,degree,density);
    end

end
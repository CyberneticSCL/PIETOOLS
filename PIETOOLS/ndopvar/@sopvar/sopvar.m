classdef (InferiorClasses={?polynomial,?dpvar})sopvar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - sopvar class
    % Represents PI maps from L_2^p[S1,S3] to L_2^q[S2,S3]
    %
    % This defines PI operators from one L2 space to another.
    %
    %   Pop: L_2^p[S1,S3] to L_2^q[S2,S3]
    %
    % Elements of this class are intended to be included in a container object.
    %
    %   The operator has the form
    %         y(S2,S3)=sum_alpha int_S1 int_S3dum  I_alpha(S_3-S_3dum) Z_d(S_2,S_3)
    %         C{\alpha} Z_d(S_3dum,S_1) x(S_3dum,S_1)
    %                         alpha \in \{0,1,-1\}^{n_3}
    %
    %   where S_1 is the variables in x(S1,S3) (input) not in S3 (common)
    %         S_2 is the variables in y(S2,S3) (output) not in S3 (common)
    %         S_3 is the variables common to x (input) and y (output)
    %         S_3dum are dummy versions of the variables in S_3
    %
    % CLASS properties
    % - P.C:   3x3x3x3...x3 Cell array of Quadpoly objects, where n_3 is the number of variables in S_3;
    %           Each quadpoly has the structure (I X Z_d(S_2,S_3)) C{alpha} (I X Z_d(S_3dum,S_1))
    %                Sachin: It is not (3^N x 1) cell, it is 3x3x...x3 where
    %                        some of the 3's are replaced by 1's when that particular
    %                        variable is not present in both domain and range.
    % - P.dims:   2x1 vector specifying the rows and columns of the operator i.e. [p,q];
    % - P.vars_S1: n1x1 cell array specifying the unique names of the spatial
    % variables in S1. There should be no duplication with S2 or S3
    % - P.vars_S2: n2x1 cell array specifying the unique names of the spatial
    % variables which appear in S2, these should be no duplication with S1 and S3
    % - P.vars_S3: n3x1 cell array specifying the unique names of the spatial
    % variables which appear in S3, these should be no duplication with S1 and
    % S2 and the naming of the dummy variables in S3dum, should be copies of
    % the variables in S3, with the appendix '_dum'
    % - P.dom_1: ordered n1x2 array specifying the spatial domain of each
    % variable in S1, so that dom_1(i,:) = [ai,bi] -> S1_i \in [ai,bi]
    % - P.dom_2: ordered n2x2 array specifying the spatial domain of each
    % variable in S2, so that dom_2(i,:) = [ai,bi] -> S2_i \in [ai,bi]
    % - P.dom_3: ordered n3x2 array specifying the spatial domain of each
    % variable in S3, so that dom_3(i,:) = [ai,bi] -> S3_i \in [ai,bi] and S3_dum_i \in [ai,bi]
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
        vars_S1 = {};   % these variables are integrated out (full integral only)
        vars_S2 = {};   % these variables are introduced (multiplier only)
        vars_S3 = {};   % these variables pass through (3 PI operators)
        dom_1 = zeros(0, 2);   % domain of s1
        dom_2 = zeros(0, 2);   % domain of s2
        dom_3 = zeros(0, 2);   % domain of s3
        dims = [1,1];          % dimension of function spaces
        params = {quadPoly(0,[],[],[1,1],{},{})};  % parameters
        % each value in the cell is a quadPoly
        % for maps to and from \R spaces, add an extra dimension?
    end

    methods
        function P = sopvar(vars3,dom3,dims,params,vars1,dom_1,vars2,dom_2)
            if nargin<=2                                                    % error handling
                error('insufficient information provided to form sopvar')
            end
            if nargin==5
                error('insufficient information provided on dom of S1')
            end
            if nargin==7
                error('insufficient information provided on dom of S2')
            end
            P.vars_S3 = vars3(:).';
            P.dims = dims;
            P.dom_3=dom3;
            if nargin<=6                                                    % S2 is empty
                P.vars_S2 = {};
                P.dom_2 = zeros(0, 2);
            end
            if nargin<=4                                                    % S1 is also empty
                P.vars_S1 = {};
                P.dom_1 = zeros(0, 2);
            end
            if nargin==3                                                    % no params specified -- create a 0-valued operator with specified domains and no S1,S2
                n= numel(vars2);                                             % total number of unique spatial variables
                if n == 0
                    P.params = {zeros(dims(1),dims(2))};
                else
                    sz=3*ones(1,n);
                    P.params = repmat({zeros(dims(1),dims(2))},sz);
                end
            end
            P.params = params;                                          % parameters specified
            if nargin>=6                                                    % parameters specified, copy S1
                P.vars_S1 = vars1(:).';
                P.dom_1=dom_1;
            end
            if nargin==8                                                % parameters specified, but no S2
                P.vars_S2 = vars2(:).';
                P.dom_2=dom_2;

            end
            % NOTE: need error checking to ensure quadpolys in params have
            % the correct format!!!!
        end
    end

    methods(Static)
        out = randsopvar(vars_S1,vars_S2,vars_S3,dim,degree,density);
    end

end
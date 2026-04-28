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
        vars = struct('in',{},'out',{});
        dom = struct('in',zeros(0,2),'out',zeros(0,2));
        dims = [1,1];
        ZL = {};
        ZR = {};
        params = {0};
    end
    properties(Access=protected)
        vars_S1;   % these variables are integrated out (full integral only)
        vars_S2;   % these variables are introduced (multiplier only)
        vars_S3;   % these variables pass through (3 PI operators)
        dom_1;   % domain of s1
        dom_2;   % domain of s2
        dom_3;   % domain of s3
    end

    methods
        function P = sopvar(params,vars,ZR,ZL,dom,dims)
            P.params = params;
            P.vars = vars;
            P.dom = dom;
            P.dims = dims;
            P.ZR = ZR; 
            P.ZL = ZL;

            P.vars_S1 = setdiff(vars.in,vars.out);
            P.vars_S2 = setdiff(vars.out,vars.in);
            P.vars_S3 = intersect(vars.in,vars.out);
            [~,idx] = ismember(P.vars_S1,vars.in);
            P.dom_1 = dom.in(idx,:);
            [~,idx] = ismember(P.vars_S2,vars.out);
            P.dom_2 = dom.out(idx,:);
            [~,idx] = ismember(P.vars_S3,vars.in);
            P.dom_3 = dom.in(idx,:);
        end
    end
    methods(Static)
        out = randsopvar(vars_S1,vars_S2,vars_S3,dim,degree,density);
        [Minv, R1Minv, R2Minv, info] = inv_1D(M, K1, K2, interval, opts);
    end

end
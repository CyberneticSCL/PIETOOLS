classdef (InferiorClasses={?opvar2d,?dpvar,?polynomial}) dopvar2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This defines the class of decision variables 
%   P: R^m0 x L2^mx x L2^my x L2^m2 to R^n0 x L2^nx x L2^ny x L2^n2
% Date: 07/12/21
% Version: 1.0
% 
% CLASS properties
% P.dim: a 4x2 array with entries [n0,m0]
%                                 [nx,mx]
%                                 [ny,my]
%                                 [n2,m2]
%
% P.R00: a n0 x m0 matrix
% P.R0x: a n0 x mx matrix valued dpvar in s1
% P.R0y: a n0 x my matrix valued dpvar in s2
% P.R02: a n0 x m2 matrix valued dpvar in s1, s2
%
% P.Rx0: a nx x m0 matrix valued dpvar in s1
% P.Rxx: 3x1 cell of nx x mx matrix valued dpvars in s1, s1_dum
% P.Rxy: a nx x my matrix valued dpvar in s1, s2
% P.Rx2: 3x1 cell of nx x m2 matrix valued dpvars in s1, s1_dum, s2_dum
%
% P.Ry0: a ny x m0 matrix valued dpvar in s2
% P.Ryx: a ny x mx matrix valued dpvar in s1, s2
% P.Ryy: 1x3 cell of ny x my matrix valued dpvars in s2, s2_dum
% P.Ry2: 1x3 cell ny x m2 matrix valued dpvars in s1, s2, s2_dum
%
% P.R20: a n2 x m0 matrix valued dpvars in s1, s2
% P.R2x: 3x1 cell of n2 x mx matrix valued dpvars in s1, s1_dum, s2
% P.R2y: 1x3 cell fo n2 x my matrix valued dpvars in s1, s2, s2_dum
% P.R22: 3x3 cell of n2 x m2 matrix valued dpvars in s1, s1_dum, s2, s2_dum
%
% P.I:    2x2 array specifying the domain as [a1, b1; a2, b2]   
% P.var1: primary polynomial variables s1 in [a1,b1], s2 in [a2,b2]
% P.var2: dummy polynomial variables s1_dum in [a1,b1], s2_dum in [a2,b2]
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - dopvar2d
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 31_01_2021  
% 06/14/2022, DJ: Update with new dim function.
% 11/30/2024, DJ: Replace default variables
%                   (ss1,ss2,tt1,tt2) --> (s1,s2,s1_dum,s2_dum);


properties
    R00 = dpvar([]);
    R0x = dpvar([]);
    R0y = dpvar([]);
    R02 = dpvar([]);
    
    Rx0 = dpvar([]);
    Rxx = {dpvar([]); dpvar([]); dpvar([])};
    Rxy = dpvar([]);
    Rx2 = {dpvar([]); dpvar([]); dpvar([])};
    
    Ry0 = dpvar([]);
    Ryx = dpvar([]);
    Ryy = {dpvar([]), dpvar([]), dpvar([])};
    Ry2 = {dpvar([]), dpvar([]), dpvar([])};
    
    R20 = dpvar([]);
    R2x = {dpvar([]); dpvar([]); dpvar([])};
    R2y = {dpvar([]), dpvar([]), dpvar([])};
    R22 = {dpvar([]), dpvar([]), dpvar([]);
           dpvar([]), dpvar([]), dpvar([]);
           dpvar([]), dpvar([]), dpvar([])};
    
    I = [0,1;0,1]
    var1 = [pvar('s1'); pvar('s2')];
    var2 = [pvar('s1_dum'); pvar('s2_dum')];
    dim = zeros(4,2);
end

properties (Dependent)
    dimdependent;
end

methods
    function [P] = dopvar2d(varargin) %constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin==1
            if ischar(varargin{1})
                if nargout==0
                    assignin('caller', varargin{1}, dopvar2d());
                end
            elseif isa(varargin{1},'double') && nargout==1
                if ~any(size(varargin{1})==[4,2])
                    error('dimension of dopvar2d must be a 4x2 integer array')
                else
                    pvar s1 s2 s1_dum s2_dum;
                    Pdim = varargin{1};     Pdom = [0,1;0,1];
                    var1 = [s1;s2];         var2 = [s1_dum;s2_dum];
                    P = dopvar2d([],Pdim,Pdom,var1,var2);
                end
            elseif isa(varargin{1},'opvar2d') && nargout==1
                P = opvar2dopvar2d(varargin{1});
                % Check if the produced object is valid
                [logval,msg] = isvalid(P);
                if ~logval
                    error(['The input is not a valid dopvar2d object: ',msg]);
                end
            elseif isa(varargin{1},'dopvar2d') && nargout==1
                P = varargin{1};
                % Check if the produced object is valid
                [logval,msg] = isvalid(P);
                if ~logval
                    error(['The input is not a valid dopvar2d object: ',msg]);
                end
            else
                error("Single input must be string or opvar2d object");
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif nargin==2
            if isa(varargin{1},'opvar2d') && (iscellstr(varargin{2}) || ischar(varargin{2}))
                P = opvar2dopvar2d(varargin{1},varargin{2});
                % Check if the produced object is valid
                [logval,msg] = isvalid(P);
                if ~logval
                    error(['A valid dopvar2d object cannot be constructed from the inputs: ',msg]);
                end
            elseif (isa(varargin{1},'double') && all(size(varargin{1})==[4,2])) && ...
                    (isa(varargin{2},'double') && all(size(varargin{2})==[2,2]))
                % Build empty dopvar2d of dimension varargin{1} on
                % domain varargin{2}
                pvar s1 s2 s1_dum s2_dum;
                Pdim = varargin{1};     Pdom = varargin{2};
                var1 = [s1;s2];         var2 = [s1_dum;s2_dum];
                P = dopvar2d([],Pdim,Pdom,var1,var2);
            elseif (isa(varargin{1},'dopvar2d') || isa(varargin{1},'opvar2d')) && isnumeric(varargin{2})
                if ~all(size(varargin{2})==[4,2])
                    error('Dimension of dopvar2d must be a 4x2 integer array')
                elseif ~all(varargin{1}.dim==varargin{2})
                    error('The size of the dopvar2d input does not match the proposed dimensions')
                else
                    if isa(varargin{1},'dopvar2d')
                        % Return the input object
                        P = varargin{1};
                        % Check if the produced object is valid
                        [logval,msg] = isvalid(P);
                        if ~logval
                            error(['The input is not a valid opvar2d object: ',msg]);
                        end
                    else
                        % Convert to dopvar
                        P = opvar2dopvar2d(varargin{1});
                    end
                end
            elseif (isa(varargin{1},'double') || isa(varargin{1},'polynomial') || isa(varargin{1},'dpvar')) && isnumeric(varargin{2})
                % Build dopvar2d from matrix varargin{1} based on
                % dimensions varargin{2}
                if ~all(size(varargin{2})==[4,2])
                    error('Dimension of opvar2d must be a 4x2 integer array')
                elseif ~(all(sum(varargin{2},1)==size(varargin{1})) || isempty(varargin{1}))
                    error('Dimension of desired opvar2d should match dimension of the input matrix')
                else
                    pvar s1 s2 s1_dum s2_dum;
                    Pdim = varargin{2};     Pdom = [0,1;0,1];
                    var1 = [s1;s2];         var2 = [s1_dum;s2_dum];
                    P = dopvar2d(varargin{1},Pdim,Pdom,var1,var2);
                end
            else
                for i=1:nargin
                    if ischar(varargin{i})
                        if nargout==0
                            assignin('caller', varargin{i}, dopvar2d());
                        end
                    else
                        error("Input must be strings");
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif nargin==3
            if (isa(varargin{1},'dopvar2d') || isa(varargin{1},'opvar2d')) && isnumeric(varargin{2})
                if ~all(size(varargin{2})==[4,2])
                    error('Dimension of dopvar2d must be a 4x2 integer array')
                elseif ~all(varargin{1}.dim==varargin{2})
                    error('The size of the dopvar2d input does not match the proposed dimensions')
                elseif ~all(varargin{1}.I==varargin{3})
                    error('The domain of the dopvar2d input does not match the proposed domain')
                else
                    if isa(varargin{1},'dopvar2d')
                        % Return the input object
                        P = varargin{1};
                        % Check if the produced object is valid
                        [logval,msg] = isvalid(P);
                        if ~logval
                            error(['The input is not a valid opvar2d object: ',msg]);
                        end
                    else
                        % Convert to dopvar
                        P = opvar2dopvar2d(varargin{1});
                    end
                end
            elseif (isa(varargin{1},'double') && all(size(varargin{1})==[4,2])) && ...
                    (isa(varargin{2},'double') && all(size(varargin{2})==[2,2])) && ...
                    (ispvar(varargin{3}) && ((prod(size(varargin{3}))==4) || (prod(size(varargin{3}))==2)))
                % Build empty dopvar2d of dimension varargin{1} on
                % domain varargin{2} in variables varargin{3}
                pvar s1_dum s2_dum;
                Pdim = varargin{1};     Pdom = varargin{2};
                var1 = [varargin{3}(1);varargin{3}(2)];
                if prod(size(varargin{3}))==4
                    var2 = [varargin{3}(3);varargin{3}(4)];
                else
                    var2 = [s1_dum;s2_dum];
                end
                P = dopvar2d([],Pdim,Pdom,var1,var2);
            elseif (isa(varargin{1},'double') || isa(varargin{1},'polynomial') || isa(varargin{1},'dpvar')) && isnumeric(varargin{2})
                % Build opvar2d from matrix varargin{1} based on
                % dimensions varargin{2}, and with domain varargin{3}.
                if ~all(size(varargin{2})==[4,2])
                    error('Dimension of opvar2d must be a 4x2 integer array')
                elseif ~(all(sum(varargin{2},1)==size(varargin{1})) || isempty(varargin{1}))
                    error('Dimension of desired opvar2d should match dimension of the input matrix')
                elseif ~isnumeric(varargin{3}) || ~all(size(varargin{3})==[2,2])
                    error('Domain should be specified as a 2x2 array')
                else
                    pvar s1 s2 s1_dum s2_dum;
                    Pdim = varargin{2};     Pdom = varargin{3};
                    var1 = [s1;s2];         var2 = [s1_dum;s2_dum];
                    P = dopvar2d(varargin{1},Pdim,Pdom,var1,var2);
                end
            else
                for i=1:nargin
                    if ischar(varargin{i})
                        if nargout==0
                            assignin('caller', varargin{i}, dopvar2d());
                        end
                    else
                        error("Input must be strings");
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif nargin==4
            if (isa(varargin{1},'dopvar2d') || isa(varargin{1},'opvar2d')) && isnumeric(varargin{2})
                if ~all(size(varargin{2})==[4,2])
                    error('Dimension of dopvar2d must be a 4x2 integer array')
                elseif ~all(varargin{1}.dim==varargin{2})
                    error('The size of the dopvar2d input does not match the proposed dimensions')
                elseif ~all(varargin{1}.I==varargin{3})
                    error('The domain of the dopvar2d input does not match the proposed domain')
                elseif ~all(size(varargin{4})==[2,2]) && ~all(size(varargin{4})==[2,1])
                    error('Variables of the dopvar2d should be specified as 2x2 pvar class object')
                elseif all(size(varargin{4})==[2,2]) && ~all(all(isequal([varargin{1}.var1,varargin{1}.var2],varargin{4})))
                    error('The variables of the dopvar2d input do not match the proposed variables')
                elseif all(size(varargin{4})==[2,1]) && ~all(isequal(varargin{1}.var1,varargin{4}))
                    error('The variables of the dopvar2d input do not match the proposed variables')
                else
                    if isa(varargin{1},'dopvar2d')
                        % Return the input object
                        P = varargin{1};
                        % Check if the produced object is valid
                        [logval,msg] = isvalid(P);
                        if ~logval
                            error(['The input is not a valid opvar2d object: ',msg]);
                        end
                    else
                        % Convert to dopvar
                        P = opvar2dopvar2d(varargin{1});
                    end
                end
            elseif (isa(varargin{1},'double') && all(size(varargin{1})==[4,2])) && ...
                    (isa(varargin{2},'double') && all(size(varargin{2})==[2,2])) && ...
                    (ispvar(varargin{3}) && (prod(size(varargin{3}))==2)) && ...
                    (ispvar(varargin{4}) && (prod(size(varargin{4}))==2))
                % Build empty dopvar2d of dimension varargin{1} on
                % domain varargin{2} in variables varargin{3}
                Pdim = varargin{1}; Pdom = varargin{2};
                var1 = [varargin{3}(1);varargin{3}(2)];
                var2 = [varargin{4}(1);varargin{4}(2)];
                P = dopvar2d([],Pdim,Pdom,var1,var2);
            elseif (isa(varargin{1},'double') || isa(varargin{1},'polynomial') || isa(varargin{1},'dpvar')) && isnumeric(varargin{2})
                % Build opvar2d from matrix varargin{1} based on
                % dimensions varargin{2}, and with domain varargin{3}.
                if ~all(size(varargin{2})==[4,2])
                    error('Dimension of opvar2d must be a 4x2 integer array')
                elseif ~(all(sum(varargin{2},1)==size(varargin{1})) || isempty(varargin{1}))
                    error('Dimension of desired opvar2d should match dimension of the input matrix')
                elseif ~isnumeric(varargin{3}) || ~all(size(varargin{3})==[2,2])
                    error('Domain should be specified as a 2x2 array')
                elseif ~ispvar(varargin{4}) || (~(prod(size(varargin{4}))==4) || ~(prod(size(varargin{4}))==2))
                    error('Variables should be specified as 2x2 pvar class object')
                else
                    pvar s1_dum s2_dum;
                    Pdim = varargin{2};     Pdom = varargin{3};
                    var1 = [varargin{4}(1);varargin{4}(2)];
                    if prod(size(varargin{4}))==4
                        var2 = [varargin{4}(3);varargin{4}(4)];
                    else
                        var2 = [s1_dum;s2_dum];
                    end
                    P = dopvar2d(varargin{1},Pdim,Pdom,var1,var2);
                end
            else
                for i=1:nargin
                    if ischar(varargin{i})
                        if nargout==0
                            assignin('caller', varargin{i}, dopvar2d());
                        end
                    else
                        error("Input must be strings");
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif nargin==5
            if (isa(varargin{1},'dopvar2d') || isa(varargin{1},'opvar2d')) && isnumeric(varargin{2})
                if ~all(size(varargin{2})==[4,2])
                    error('Dimension of dopvar2d must be a 4x2 integer array')
                elseif ~all(varargin{1}.dim==varargin{2})
                    error('The size of the dopvar2d input does not match the proposed dimensions')
                elseif ~all(varargin{1}.I==varargin{3})
                    error('The domain of the dopvar2d input does not match the proposed domain')
                elseif ~all(size(varargin{4})==[2,1]) || ~all(size(varargin{5})==[2,1])
                    error('Variables of the dopvar2d should be specified as two 2x1 pvar class object')
                elseif ~all(isequal(varargin{1}.var1,varargin{4}))
                    error('The primary variables of the dopvar2d input do not match the proposed variables')
                elseif ~all(isequal(varargin{1}.var2,varargin{5}))
                    error('The secondary variables of the dopvar2d input do not match the proposed variables')
                else
                    if isa(varargin{1},'dopvar2d')
                        % Return the input object
                        P = varargin{1};
                        % Check if the produced object is valid
                        [logval,msg] = isvalid(P);
                        if ~logval
                            error(['The input is not a valid opvar2d object: ',msg]);
                        end
                    else
                        % Convert to dopvar
                        P = opvar2dopvar2d(varargin{1});
                    end
                end
            elseif (isa(varargin{1},'double') || isa(varargin{1},'polynomial') || isa(varargin{1},'dpvar')) && isnumeric(varargin{2})
                % Build dopvar2d from matrix varargin{1} based on
                % dimensions varargin{2}, and with domain varargin{3}.
                if ~all(size(varargin{2})==[4,2])
                    error('Dimension of opvar2d must be a 4x2 integer array')
                elseif ~(all(sum(varargin{2},1)==size(varargin{1})) || isempty(varargin{1}))
                    error('Dimension of desired opvar2d should match dimension of the input matrix')
                elseif ~isnumeric(varargin{3}) || ~all(size(varargin{3})==[2,2])
                    error('Domain should be specified as a 2x2 array')
                elseif ~ispvar(varargin{4}) || ~(prod(size(varargin{4}))==2) || ...
                        ~ispvar(varargin{5}) || ~(prod(size(varargin{5}))==2)
                    error('Variables should be specified as two 2x1 pvar class object')
                else
                    Pdim = varargin{2};
                    P = dopvar2d();
                    P.dim = Pdim;
                    P = set(P,'I',varargin{3});
                    P = set(P,'var1',[varargin{4}(1);varargin{4}(2)],'nosubs');
                    P = set(P,'var2',[varargin{5}(1);varargin{5}(2)],'nosubs');
                    
                    % Set the parameters based on value of first argument
                    if ~isempty(varargin{1})
                        mat = varargin{1};
                        rind = cumsum(P.dim(:,1));
                        cind = cumsum(P.dim(:,2));
                        
                        P = set(P,'R00',dpvar(mat(1:rind(1),1:cind(1))));
                        P = set(P,'Rx0',dpvar(mat(rind(1)+1:rind(2),1:cind(1))));
                        P = set(P,'Ry0',dpvar(mat(rind(2)+1:rind(3),1:cind(1))));
                        P = set(P,'R20',dpvar(mat(rind(3)+1:rind(4),1:cind(1))));
                        
                        P = set(P,'R0x',dpvar(mat(1:rind(1),cind(1)+1:cind(2))));
                        L(1).type = '.';    L(2).type = '{}';
                        L(1).subs = 'Rxx';  L(2).subs = {1};
                        P = subsasgn(P,L,dpvar(mat(rind(1)+1:rind(2),cind(1)+1:cind(2))));
                        P = set(P,'Ryx',dpvar(mat(rind(2)+1:rind(3),cind(1)+1:cind(2))));
                        L(1).subs = 'Rx2';
                        P = subsasgn(P,L,dpvar(mat(rind(3)+1:rind(4),cind(1)+1:cind(2))));
                        
                        P = set(P,'R0y',dpvar(mat(1:rind(1),cind(2)+1:cind(3))));
                        P = set(P,'Rxy',dpvar(mat(rind(1)+1:rind(2),cind(2)+1:cind(3))));
                        L(1).subs = 'Ryy';
                        P = subsasgn(P,L,dpvar(mat(rind(2)+1:rind(3),cind(2)+1:cind(3))));
                        L(1).subs = 'Ry2';
                        P = subsasgn(P,L,dpvar(mat(rind(3)+1:rind(4),cind(2)+1:cind(3))));
                        
                        set(P,'R02',dpvar(mat(1:rind(1),cind(3)+1:cind(4))));
                        L(1).subs = 'R2x';
                        P = subsasgn(P,L,dpvar(mat(rind(1)+1:rind(2),cind(3)+1:cind(4))));
                        L(1).subs = 'R2y';
                        P = subsasgn(P,L,dpvar(mat(rind(2)+1:rind(3),cind(3)+1:cind(4))));
                        L(1).subs = 'R22';
                        P = subsasgn(P,L,dpvar(mat(rind(3)+1:rind(4),cind(3)+1:cind(4))));
                        
                        % Check if the produced object is valid
                        [logval,msg] = isvalid(P);
                        if ~logval
                            %fprintf(2,['\n > > >  Warning: ',msg,'  < < < \n'])
                            error(['A valid opvar2d object cannot be constructed from the provided inputs: ',msg]);
                        end
                    end
                end
            else
                for i=1:nargin
                    if ischar(varargin{i})
                        if nargout==0
                            assignin('caller', varargin{i}, dopvar2d());
                        end
                    else
                        error("Input must be strings");
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            for i=1:nargin
                if ischar(varargin{i})
                    if nargout==0
                        assignin('caller', varargin{i}, dopvar2d());
                    end
                else
                    error("Input must be strings");
                end
            end
        end
    end
    % % % =========================================================== % % %
    
    function [obj] = set.I(obj,val)
        if ~all(size(val)==[2,2]) || ~all(val(:,1)<val(:,2))
            error('The opvar2d domain "I" should be specified as a 2x2 array, with values in the first column strictly smaller than those in the second')
        else
            obj.I = val;
        end
    end
    
    function [obj] = set.dim(obj,Pdim_new)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % obj = set.dim(obj,Pdim_new) takes in one opvar2d object Pop, and returns
    % a dopvar2d object out with out.dim = Pdim_new, if possible.
    %
    % Version 1.0
    % Date: 06/13/22
    %
    % INPUT
    % Pop:  dopvar2d class object;
    % dim:  a 4x2 array of nonnegative integers, specifying desired dimensions
    %       as [n0, m0; nx, mx; ny, my; n2, m2].
    %
    % OUTPUT
    % out:  dopvar2d object with Pop.dim = Pdim_new. so that
    %       Pop: R^m0 x L2^mx[x] x L2^my[y] x L2^m2[x,y]
    %                           --> R^n0 x L2^nx[x] x L2^ny[y] x L2^n2[x,y].
    %       Note: If Pop.dim already has the desired dimensions, then out=Pop.
    %       Otherwise, the altered dimensions must correspond to parameters
    %       which are all zeroes. These parameters will be adjusted to match
    %       the new desired dimensions, so that e.g. [0,0;0,0] may become
    %       [0;0;0]. Any nonzero parameters of Pop will not be adjusted, and if
    %       their dimension does not match the desired dimensions, an error is
    %       produced.
    %
    % NOTES:
    % For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
    % or D. Jagt at djagt@asu.edu
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % If you modify this code, document all changes carefully and include date
    % authorship, and a brief description of modifications
    %
    % Initial coding DJ - 13/06/2022
        
        % Make sure the specified dimensions makes sense
        if ~all(size(Pdim_new)==[4,2]) || any(Pdim_new(:)<0)
            error('Dimension of an opvar2d object must be specified as a 4x2 array of nonnegative integers')
        end
        
        % If we're initializing a dopvar2d object, don't bother checking
        % parameter sizes, but simply set the new parameters as zeros.
        Pdim = obj.dim;
        if all(Pdim(:)==0)
            Rparams = {'R00', 'R0x', 'R0y', 'R02';
                       'Rx0', 'Rxx', 'Rxy', 'Rx2';
                       'Ry0', 'Ryx', 'Ryy', 'Ry2';
                       'R20', 'R2x', 'R2y', 'R22'};
            for k=1:numel(Rparams)
                PR = obj.(Rparams{k});
                [rnum,cnum] = ind2sub(size(Rparams),k);
                PRdim_new = [Pdim_new(rnum,1),Pdim_new(cnum,2)];
                if all(PRdim_new==0)
                    continue
                elseif ~isa(PR,'cell')
                    % Set all-zero parameter, without additional checks
                    obj = set(obj,Rparams{k},dpvar(zeros(PRdim_new)),'nocheck');
                elseif isa(PR,'cell')
                    PR_new = cell(size(PR));
                    for l=1:numel(PR)
                        PR_new{l} = dpvar(zeros(PRdim_new));
                    end
                    obj = set(obj,Rparams{k},PR_new,'nocheck');
                end
            end
            return
        end
        
        % Compare the proposed dimensions to the current dimensions
        Pdim = obj.dim;
        dim_match = Pdim_new==Pdim;
        if all(dim_match(:))    % if all dimensions match, do nothing
            return
        else
            % % If some dimensions do not match, see if these dimensions correspond
            % % to parameters which are just zero. We allow the dimensions of these
            % % zero parameters to be adjusted.
            
            % Initialize a new opvar2d object with the proposed dimensions, and
            % otherwise the same properties as the input operator.
            
            % Next, assign the parameters of this operator
            Rparams = {'R00', 'R0x', 'R0y', 'R02';
                       'Rx0', 'Rxx', 'Rxy', 'Rx2';
                       'Ry0', 'Ryx', 'Ryy', 'Ry2';
                       'R20', 'R2x', 'R2y', 'R22'};
            for k=1:numel(Rparams)
                % Extract the parameter, and establish the dimensions it should
                % have in the output operator.
                PR = obj.(Rparams{k});
                [rnum,cnum] = ind2sub(size(Rparams),k);
                PRdim_new = [Pdim_new(rnum,1),Pdim_new(cnum,2)];
                
                % If the size of the parameter does not match the proposed
                % dimensions, check if the parameter is all zeroes. If so,
                % expand/compress this parameter to the appropriate dimensions.
                if ~isa(PR,'cell') && all(size(PR)==PRdim_new)
                    % The dimensions match.
                    obj.(Rparams{k}) = obj.(Rparams{k});
                elseif isempty(PR) || ...
                        (isa(PR,'double') && all(all(PR==0))) || ...
                        (isa(PR,'polynomial') && all(all(isequal(PR,0)))) || ...
                        (isa(PR,'dpvar') && all(all(PR.C==0)))
                    % The parameter is all zeroes
                    obj.(Rparams{k}) = dpvar(zeros(PRdim_new));
                elseif isa(PR,'cell')
                    PR_new = cell(size(PR));
                    for l=1:numel(PR)
                        PRR = PR{l};
                        if all(size(PRR)==[Pdim_new(rnum,1),Pdim_new(cnum,2)])
                            % The dimensions match.
                            PR_new{l} = PR{l};
                        elseif isempty(PRR) || ...
                                (isa(PRR,'double') && all(all(PRR==0))) || ...
                                (isa(PRR,'polynomial') && all(all(isequal(PRR,0)))) || ...
                                (isa(PRR,'dpvar') && all(all(PRR.C==0)))
                            % The parameter is all zeroes
                            PR_new{l} = dpvar(zeros(PRdim_new));
                        else
                            error(['Cannot adjust the dimensions of the opvar2d object; the proposed dimensions do not match those of the parameter "',(Rparams{k}),'{',num2str(l),'}".']);
                        end
                    end
                    obj.(Rparams{k}) = PR_new;
                else
                    error(['Cannot adjust the dimensions of the opvar2d object; the proposed dimensions do not match those of the parameter "',(Rparams{k}),'".']);
                end
            end
        end
        
    end
    
    
    % % % =========================================================== % % %    
    
    function [d] = get.dimdependent(obj)
        % For consistent dimensions, the values in each of the  following 
        % vectors should be identical.
        N0 = [size(obj.R00,1); size(obj.R0x,1); size(obj.R0y,1); size(obj.R02,1)];
        Nx = [size(obj.Rx0,1); size(obj.Rxy,1);
            size(obj.Rxx{1},1); size(obj.Rxx{2},1); size(obj.Rxx{3},1);
            size(obj.Rx2{1},1); size(obj.Rx2{2},1); size(obj.Rx2{3},1)];
        Ny = [size(obj.Ry0,1); size(obj.Ryx,1);
            size(obj.Ryy{1},1); size(obj.Ryy{2},1); size(obj.Ryy{3},1);
            size(obj.Ry2{1},1); size(obj.Ry2{2},1); size(obj.Ry2{3},1)];
        N2 = [size(obj.R20,1);
            size(obj.R2x{1},1); size(obj.R2x{2},1); size(obj.R2x{3},1);
            size(obj.R2y{1},1); size(obj.R2y{2},1); size(obj.R2y{3},1);
            size(obj.R22{1,1},1); size(obj.R22{1,2},1); size(obj.R22{1,3},1);
            size(obj.R22{2,1},1); size(obj.R22{2,2},1); size(obj.R22{2,3},1);
            size(obj.R22{3,1},1); size(obj.R22{3,2},1); size(obj.R22{3,3},1)];
        
        M0 = [size(obj.R00,2); size(obj.Rx0,2); size(obj.Ry0,2); size(obj.R20,2)];
        Mx = [size(obj.R0x,2); size(obj.Ryx,2);
            size(obj.Rxx{1},2); size(obj.Rxx{2},2); size(obj.Rxx{3},2);
            size(obj.R2x{1},2); size(obj.R2x{2},2); size(obj.R2x{3},2)];
        My = [size(obj.R0y,2); size(obj.Rxy,2);
            size(obj.Ryy{1},2); size(obj.Ryy{2},2); size(obj.Ryy{3},2);
            size(obj.R2y{1},2); size(obj.R2y{2},2); size(obj.R2y{3},2)];
        M2 = [size(obj.R02,2);
            size(obj.Rx2{1},2); size(obj.Rx2{2},2); size(obj.Rx2{3},2);
            size(obj.Ry2{1},2); size(obj.Ry2{2},2); size(obj.Ry2{3},2);
            size(obj.R22{1,1},2); size(obj.R22{1,2},2); size(obj.R22{1,3},2);
            size(obj.R22{2,1},2); size(obj.R22{2,2},2); size(obj.R22{2,3},2);
            size(obj.R22{3,1},2); size(obj.R22{3,2},2); size(obj.R22{3,3},2)];
        
        % If the output dimensions of the parameters along any row match,
        % set this output dimension as the dimension of the operator.
        % Otherwise, set the dimension to NaN.
        N0 = N0(N0~=0);
        if isempty(N0)
            n0=0;
        elseif all(N0/max(N0)==1)
            n0 = max(N0);
        else
            n0 = nan;
        end
        Nx = Nx(Nx~=0);
        if isempty(Nx)
            nx=0;
        elseif all(Nx/max(Nx)==1)
            nx = max(Nx);
        else
            nx = nan;
        end
        Ny = Ny(Ny~=0);
        if isempty(Ny)
            ny=0;
        elseif all(Ny/max(Ny)==1)
            ny = max(Ny);
        else
            ny = nan;
        end
        N2 = N2(N2~=0);
        if isempty(N2)
            n2=0;
        elseif all(N2/max(N2)==1)
            n2 = max(N2);
        else
            n2 = nan;
        end
        
        % If the input dimensions of the parameters along any column match,
        % set this input dimension as the dimension of the operator.
        % Otherwise, set the dimension to NaN.
        M0 = M0(M0~=0);
        if isempty(M0)
            m0=0;
        elseif all(M0/max(M0)==1)
            m0 = max(M0);
        else
            m0 = nan;
        end
        Mx = Mx(Mx~=0);
        if isempty(Mx)
            mx=0;
        elseif all(Mx/max(Mx)==1)
            mx = max(Mx);
        else
            mx = nan;
        end
        My = My(My~=0);
        if isempty(My)
            my=0;
        elseif all(My/max(My)==1)
            my = max(My);
        else
            my = nan;
        end
        M2 = M2(M2~=0);
        if isempty(M2)
            m2=0;
        elseif all(M2/max(M2)==1)
            m2 = max(M2);
        else
            m2 = nan;
        end
        
        d = [n0 m0; nx mx; ny my; n2 m2];
    end
    
    function [val] = get.dim(obj)
        val = obj.dimdependent;
    end
end
end
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

% P.R00: a n0 x m0 matrix
% P.R0x: a n0 x mx matrix valued dpvar in ss1
% P.R0y: a n0 x my matrix valued dpvar in ss2
% P.R02: a n0 x m2 matrix valued dpvar in ss1, ss2

% P.Rx0: a nx x m0 matrix valued dpvar in ss1
% P.Rxx: 3x1 cell of nx x mx matrix valued dpvars in ss1, tt1
% P.Rxy: a nx x my matrix valued dpvar in ss1, ss2
% P.Rx2: 3x1 cell of nx x m2 matrix valued dpvars in ss1, tt1, ss2

% P.Ry0: a ny x m0 matrix valued dpvar in ss2
% P.Ryx: a ny x mx matrix valued dpvar in ss1, ss2
% P.Ryy: 1x3 cell of ny x my matrix valued dpvars in ss2, tt2
% P.Ry2: 1x3 cell ny x m2 matrix valued dpvars in ss1, ss2, tt2

% P.R20: a n2 x m0 matrix valued dpvars in ss1, ss2
% P.R2x: 3x1 cell of n2 x mx matrix valued dpvars in ss1, tt1, ss2
% P.R2y: 1x3 cell fo n2 x my matrix valued dpvars in ss1, ss2, tt2
% P.R22: 3x3 cell of n2 x m2 matrix valued dpvars in ss1, tt1, ss2, tt2

% P.I: domain    
% P.var1: polynomial variables ss1, ss2, (internal property, recommend not to modify)
% P.var2: polynomial variables tt1, tt2 (internal property, recommend not to modify)
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - opvar2d
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
%   ^ Based heavily on "@opvar"-opvar code by SS ^



    properties
        R00 = [];
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
        var1 = [pvar('ss1'); pvar('ss2')];
        var2 = [pvar('tt1'); pvar('tt2')];
        dim = zeros(4,2);
    end
    
    properties (Dependent)
        dimdependent;
    end
    
    methods
        function [P] = dopvar2d(varargin) %constructor
            if nargin==1
                if ischar(varargin{1})
                    if nargout==0
                        assignin('caller', varargin{1}, dopvar2d());
                    end
                elseif isa(varargin{1},'double') && nargout==1
                    if ~any(size(varargin{1})==[4,2])
                        error('dimension of dopvar2d must be a 4x2 integer array')
                    else
                        P = dopvar2d();
                        P.dim = varargin{1};
                    end
                elseif isa(varargin{1},'opvar2d') && nargout==1
                    P = opvar2dopvar2d(varargin{1});
                elseif isa(varargin{1},'dopvar2d') && nargout==1
                    P = varargin{1};
                else
                    error("Single input must be string or opvar2d object");
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
        end
        
        function [obj] = set.R00(obj,R00) 
            obj.R00 = R00;
        end
        function [obj] = set.R0x(obj,R0x) 
            obj.R0x = R0x;
        end
        function [obj] = set.R0y(obj,R0y) 
            obj.R0y = R0y;
        end
        function [obj] = set.R02(obj,R02) 
            obj.R02 = R02;
        end
        
        function [obj] = set.Rx0(obj,Rx0) 
            obj.Rx0 = Rx0;
        end
        function [obj] = set.Rxx(obj,Rxx) 
            obj.Rxx = Rxx;
        end
        function [obj] = set.Rxy(obj,Rxy) 
            obj.Rxy = Rxy;
        end
        function [obj] = set.Rx2(obj,Rx2) 
            obj.Rx2 = Rx2;
        end
        
        function [obj] = set.Ry0(obj,Ry0) 
            obj.Ry0 = Ry0;
        end
        function [obj] = set.Ryx(obj,Ryx) 
            obj.Ryx = Ryx;
        end
        function [obj] = set.Ryy(obj,Ryy) 
            obj.Ryy = Ryy;
        end
        function [obj] = set.Ry2(obj,Ry2) 
            obj.Ry2 = Ry2;
        end
        
        function [obj] = set.R20(obj,R20) 
            obj.R20 = R20;
        end
        function [obj] = set.R2x(obj,R2x) 
            obj.R2x = R2x;
        end
        function [obj] = set.R2y(obj,R2y) 
            obj.R2y = R2y;
        end
        function [obj] = set.R22(obj,R22) 
            obj.R22 = R22;
        end
        
        function [d] = get.dimdependent(obj)
            % for consistent dimensions following vectors should have
            % all values in each of them to be equal or zero.
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
        function [obj] = set.dim(obj,val)
            obj.dim = val;
            if isempty(obj.R00)
                obj.R00 = zeros(val(1,:));
            end
            if isempty(obj.R0x)
                obj.R0x = zeros(val(1,1),val(2,2));
            end
            if isempty(obj.R0y)
                obj.R0y = zeros(val(1,1),val(3,2));
            end
            if isempty(obj.R02)
                obj.R02 = zeros(val(1,1),val(4,2));
            end
            if isempty(obj.Rx0)
                obj.Rx0 = zeros(val(2,1),val(1,2));
            end
            if isempty(obj.Rxy)
                obj.Rxy = zeros(val(2,1),val(3,2));
            end
            if isempty(obj.Ry0)
                obj.Ry0 = zeros(val(3,1),val(1,2));
            end
            if isempty(obj.Ryx)
                obj.Ryx = zeros(val(3,1),val(2,2));
            end
            if isempty(obj.R20)
                obj.R20 = zeros(val(4,1),val(1,2));
            end
            
            for i=1:3
                if isempty(obj.Rxx{i,1})
                    obj.Rxx{i,1} = zeros(val(2,:));
                end          
                if isempty(obj.Rx2{i,1})
                    obj.Rx2{i,1} = zeros(val(2,1),val(4,2));
                end
                
                if isempty(obj.Ryy{1,i})
                    obj.Ryy{1,i} = zeros(val(3,:));
                end
                if isempty(obj.Ry2{1,i})
                    obj.Ry2{1,i} = zeros(val(3,1),val(4,2));
                end           
            
                if isempty(obj.R2x{i,1})
                    obj.R2x{i,1} = zeros(val(4,1),val(2,2));
                end
                if isempty(obj.R2y{1,i})
                    obj.R2y{1,i} = zeros(val(4,1),val(3,2));
                end
                
                for j=1:3
                    if isempty(obj.R22{i,j})
                        obj.R22{i,j} = zeros(val(4,:));
                    end
                end
            end
            
        end
    end
end
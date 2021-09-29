function [prog,Pop] = poslpivar_2d(prog,n,I,d,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop,LLL,bZ_col] = poslpivar_2d(prog,n,I,d,options) declares 
% a positive 0112D PI operator. 
%
% Pop = Z' * T * Z
% T = [ T_{00,00} T_{00,x1} T_{00,x2} T_{00,x3} T_{00,y1} T_{00,y2} T_{00,y3} T_{00,21} ... T_{00,29}]
%     [      :         :         :         :         :         :         :         :             :   ] > 0
%     [ T_{29,00} T_{29,x1} T_{29,x2} T_{29,x3} T_{29,y1} T_{29,y2} T_{29,y3} T_{29,21} ... T_{29,29}]
%
% Z*u = [                                 u0                ]
%       [                          Zxo(x)*ux(x)             ]
%       [               int_a^x Zxa(x,ss)*ux(ss) dss        ]
%       [               int_x^b Zxb(x,ss)*ux(ss) dss        ]
%       [                          Zyo(x)*uy(y)             ]
%       [               int_c^y Zya(y,tt)*uy(tt) dtt        ]
%       [               int_y^d Zyb(y,tt)*uy(tt) dtt        ]
%       [                       Z2oo(x,y)*u2(x,y)           ]
%       [            int_a^x Z2ao(x,y,ss)*u2(ss,y) dss      ]
%       [            int_x^b Z2bo(x,y,ss)*u2(ss,y) dss      ]
%       [            int_c^y Z2oa(x,y,tt)*u2(x,tt) dtt      ]
%       [            int_y^d Z2ob(x,y,tt)*u2(x,tt) dtt      ]
%       [ int_a^x int_c^y Z2aa(x,y,ss,tt)*u2(ss,tt) dtt dss ]
%       [ int_a^x int_y^d Z2ab(x,y,ss,tt)*u2(ss,tt) dtt dss ]
%       [ int_x^b int_c^y Z2ba(x,y,ss,tt)*u2(ss,tt) dtt dss ]
%       [ int_x^b int_y^d Z2bb(x,y,ss,tt)*u2(ss,tt) dtt dss ]
%
% where each function Z.. =  Z_d(x,y,tt,ss) \otimes I_n.. with Z_d(x) a 
% vector of monomials in variables x,y,tt,ss of total degree d or less. 
% 
% INPUT 
%   prog: SOS program to modify.
%   n(1): dimension of real part
%   n(2): dimension of L2[x] part
%   n(3): dimension of L2[y] part
%   n(4): dimension of L2[x,y] part

%   I = [l u] interval of integration (a=I(1,1), b=I(1,2), c=I(2,1), d=I(2,2))
%   -Optional INPUTS
%   d, structure with fields dx, dy, d2, where:
%   dx{1}: degree of s in Zx^o(s)
%   dx{2}(1): degree of s in Zx^a(s,th), defaults to d(1)
%   dx{2}(2): degree of th in Zx^a(s,th), defaults to d(1)
%   dx{2}(3): joint degree of s,th in Zx^a(s,th), defaults to d(2,1)+d(2,2)
%   dx{3}(1): degree of s in Zx^b(s,th), defaults to d(1)
%   dx{3}(2): degree of th in Zx^b(s,th), defaults to d(1)
%   dx{3}(3): joint degree of s,th in Zx^b(s,th), defaults to d(3,1)+d(3,2)

%   d2{1,1}(1): degree of s1 in Z2^oo(s1,s2)
%   d2{1,1}(2): degree of s2 in Z2^oo(s1,s2)
%   d2{1,1}(3): joint degree of s1,s2 in Z2^oo(s1,s2), defaults to d{1,1}(1)+d{1,1}(2)

%   d2{2,1}(1): degree of s1 in Z2^ao(s1,th1,s2), defaults to d2{1,1}(1)
%   d2{2,1}(2): degree of th1 in Z2^ao(s1,th1,s2), defaults to d2{1,1}(1)
%   d2{2,1}(3): joint degree of s1,th1 in Z2^ao(s1,th1,s2), defaults to d{2,1}(1)+d{2,1}(2)
%   d2{2,1}(4): degree of s2 in Z2^ao(s1,th1,s2), defaults to d2{1,1}(2)
%   d2{2,1}(5): joint degree of s1,th1,s2 in Z2^ao(s1,th1,s2), defaults to d{2,1}(3)+d{2,1}(4)
%   d2{3,1}(1): degree of s1 in Z2^bo(s1,th1,s2), defaults to d2{1,1}(1)
%   d2{3,1}(2): degree of th1 in Z2^bo(s1,th1,s2), defaults to d2{1,1}(1)
%   d2{3,1}(3): joint degree of s1,th1 in Z2^bo(s1,th1,s2), defaults to d{3,1}(1)+d{3,1}(2)
%   d2{3,1}(4): degree of s2 in Z2^bo(s1,th1,s2), defaults to d2{1,1}(2)
%   d2{3,1}(5): joint degree of s1,th1,s2 in Z2^bo(s1,th1,s2), defaults to d{3,1}(3)+d{3,1}(4)

%   d2{1,2}(1): degree of s1 in Z2^oa(s1,s2,th2), defaults to d2{1,1}(1)
%   d2{1,2}(2): degree of s2 in Z2^oa(s1,s2,th2), defaults to d2{1,1}(2)
%   d2{1,2}(3): degree of th2 in Z2^oa(s1,s2,th2), defaults to d2{1,1}(2)
%   d2{1,2}(4): joint degree of s2,th2 in Z2^oa(s1,s2,th2), defaults to d2{1,2}(2)+d{1,2}(3)
%   d2{1,3}(1): degree of s1 in Z2^ob(s1,s2,th2), defaults to d2{1,1}(1)
%   d2{1,3}(2): degree of s2 in Z2^ob(s1,s2,th2), defaults to d2{1,1}(2)
%   d2{1,3}(3): degree of th2 in Z2^ob(s1,s2,th2), defaults to d2{1,1}(2)
%   d2{1,3}(4): joint degree of s2,th2 in Z2^ob(s1,s2,th2), defaults to d2{1,3}(2)+d{1,3}(3)

%   d2{2,2}(1,1): degree of s1 in Z2^aa(s1,th1,s2,th2), defaults to d2{1,1}(1)
%   d2{2,2}(2,1): degree of th1 in Z2^aa(s1,th1,s2,th2), defaults to d2{1,1}(1)
%   d2{2,2}(3,1): joint degree of s1,th1 in Z2^aa(s1,th1,s2,th2), defaults to d{2,2}(1,1)+d{2,2}(2,1)
%   d2{2,2}(1,2): degree of s2 in Z2^aa(s1,th1,s2,th2), defaults to d2{1,1}(2)
%   d2{2,2}(2,2): degree of th2 in Z2^aa(s1,th1,s2,th2), defaults to d2{1,1}(2)
%   d2{2,2}(3,2): joint degree of s2,th2 in Z2^aa(s1,th1,s2,th2), defaults to d{2,2}(1,2)+d{2,2}(2,2)

%   options.psatz=1 if this is a psatz term. options.psatz=0 otherwise
%   options.exclude is a length 16 binary vector where 
%      options.exclude(i)=1 if we want to set $T_{ij}=0$ for j=1...16
%   options.sep is a length 5 binary vector where
%      options.sep(1) = 1 if Rxx{2} = Rxx{3}
%      options.sep(2) = 1 if Ryy{2} = Ryy{3}
%      options.sep(3) = 1 if R22{2,1} = R22{3,1}
%      options.sep(4) = 1 if R22{1,2} = R22{1,3}
%      options.sep(5) = 1 if R22{2,2} = R22{3,2} = R22{2,3} = R22{3,3}
% 
% OUTPUT 
%   prog: modified SOS program with new variables and constraints
%   Pop: operator structure
%   - OPTIONAL OUTPUT - CURRENTLY NOT SUPPORTED!!!
%   LLL: positive matrix variable used to create Pop
%   bZ_col: Monomial function used to create Pop
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - poslpivar_2d
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
% Initial coding DJ - 09_28_2021

% % % Set-up % % %

% % Extract the input arguments
switch nargin
    case 2
        error('Not enough inputs!')
    case 3
        dx = {1;[1;1;1];[1;1;1]};
        dy = {1,[1;1;1],[1;1;1]};
        d2 = {[1;1;2],[1;1;1;1;2],[1;1;1;1;2];
              [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1];
              [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1]};
          options.psatz = 0;
          options.exclude = zeros(1,16);
          options.diag = 0;
          options.sep = zeros(1,6);
    case 4
        if ~isfield(d,'dx') && ~isfield(d,'dy') && ~isfield(d,'d2')
            fprintf('\n Warning: No degrees are specified. Continuing with default values. \n')
            dx = {1;[1;1;1];[1;1;1]};
            dy = {1,[1;1;1],[1;1;1]};
            d2 = {[1;1;2],[1;1;1;1;2],[1;1;1;1;2];
                [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1];
                [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1]};
        else
            if isfield(d,'dx')
                dx = d.dx;
            elseif isfield(d,'dy')
                dy = d.dy;
                dx = dy;
            else
                d2 = d.d2;
                dx = cell(3,1);
                dx{1} = d2{1,1}(1);  dx{2} = d2{2,1}(1:3);  dx{3} = d2{3,1}(1:3);
            end
            if isfield(d,'dy')
                dy = d.dy;
            elseif isfield(d,'dx')
                dy = dx;
            else
                d2 = d.d2;
                dy = cell(1,3);
                dy{1} = d2{1,1}(2);  dy{2} = d2{1,2}(2:4);  dy{3} = d2{1,3}(2:4);
            end
            if isfield(d,'d2')
                d2 = d.d2;
            else
                d2 = {[dx{1};dy{1}],[dx{1};dy{2}(:);ceil(0.5*(dx{1}+dy{2}(end)))],[dx{1};dy{3}(:);ceil(0.5*(dx{1}+dy{3}(end)))];
                    [dx{2}(:);dy{1};ceil(0.5*(dx{2}(end)+dy{1}))],[dx{2};dy{2};ceil(0.5*(dx{2}+dy{2}))],[dx{2};dy{3};ceil(0.5*(dx{2}+dy{3}))];
                    [dx{3}(:);dy{1};ceil(0.5*(dx{2}(end)+dy{1}))],[dx{3};dy{2};ceil(0.5*(dx{3}+dy{2}))],[dx{3};dy{3};ceil(0.5*(dx{3}+dy{3}))]};
            end
        end
        options.psatz = 0;
        options.exclude = zeros(1,16);
        options.diag = 0;
        options.sep = zeros(1,6);
    case 5
        if ~isfield(d,'dx') && ~isfield(d,'dy') && ~isfield(d,'d2')
            fprintf('\n Warning: No degrees are specified. Continuing with default values. \n')
            dx = {1;[1;1;1];[1;1;1]};
            dy = {1,[1;1;1],[1;1;1]};
            d2 = {[1;1;2],[1;1;1;1;2],[1;1;1;1;2];
                [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1];
                [1;1;1;1;2],[1,1,1;1,1,1;1,1,1],[1,1,1;1,1,1;1,1,1]};
        else
            if isfield(d,'dx')
                dx = d.dx;
            elseif isfield(d,'dy')
                dy = d.dy;
                dx = dy;
            else
                d2 = d.d2;
                dx = cell(3,1);
                dx{1} = d2{1,1}(1);  dx{2} = d2{2,1}(1:3);  dx{3} = d2{3,1}(1:3);
            end
            if isfield(d,'dy')
                dy = d.dy;
            elseif isfield(d,'dx')
                dy = dx;
            else
                d2 = d.d2;
                dy = cell(1,3);
                dy{1} = d2{1,1}(2);  dy{2} = d2{1,2}(2:4);  dy{3} = d2{1,3}(2:4);
            end
            if isfield(d,'d2')
                d2 = d.d2;
            else
                d2 = {[dx{1};dy{1}],[dx{1};dy{2}(:);ceil(0.5*(dx{1}+dy{2}(end)))],[dx{1};dy{3}(:);ceil(0.5*(dx{1}+dy{3}(end)))];
                    [dx{2}(:);dy{1};ceil(0.5*(dx{2}(end)+dy{1}))],[dx{2};dy{2};ceil(0.5*(dx{2}+dy{2}))],[dx{2};dy{3};ceil(0.5*(dx{2}+dy{3}))];
                    [dx{3}(:);dy{1};ceil(0.5*(dx{2}(end)+dy{1}))],[dx{3};dy{2};ceil(0.5*(dx{3}+dy{2}))],[dx{3};dy{3};ceil(0.5*(dx{3}+dy{3}))]};
            end
        end
        if ~isfield(options,'psatz')
            options.psatz=0;
        end
        if ~isfield(options,'exclude')
            options.exclude = zeros(1,16);
        end
        if isfield(options,'diag') && options.diag==1
            fprintf(2,'Warning: ''diag'' option is not supported for 2D PDEs, ignoring this input.'); 
        end
        options.diag=0;
        if ~isfield(options,'sep')
            options.sep = zeros(1,6);
        end
end
if length(n)==1
    warning('Only 1 dimension is provided, assuming this to refer to the PDE state')
    n = [0,0,0,n(2)];
elseif length(n)==2
    warning('Only 2 dimensions are provided, assuming these to refer to the ODE and PDE states')
    n = [n(1),0,0,n(2)];
elseif length(n)~=4
    error('n must be a length 4 vector')
end

if any(I(:,1)>=I(:,2))
    error('I(:,1) must be less than I(:,2)')
end

% Specify the degrees of the monomials
if ~iscell(dx)
    error('dx is a 3-cell structure')
end
if ~iscell(dy)
    error('dy is a 3-cell structure')
end
if ~iscell(d2)
    error('d2 is a 3x3-cell structure')
end

if length(dx(:))==1
    dx{2}=[dx{1},dx{1},2*dx{1}];
    dx{3}=dx{2};
elseif length(dx(:))==2
    if length(dx{2})==1
        dx{2}(2)=dx{2}(1);
        dx{2}(3)=dx{2}(1);
    elseif length(dx{2})==2
        dx{2}(3) = ceil(0.5*(dx{2}(1) + dx{2}(2)));
    end
    dx{3}=dx{2};
else
    if length(dx{2})==1
        dx{2}(2)=dx{2}(1);
        dx{2}(3)=dx{2}(1);
    elseif length(dx{2})==2
        dx{2}(3) = ceil(0.5*(dx{2}(1) + dx{2}(2)));
    end
    if length(dx{3})==1
        dx{3}(2)=dx{3}(1);
        dx{3}(3)=dx{3}(1);
    elseif length(dx{3})==2
        dx{3}(3) = ceil(0.5*(dx{3}(1) + dx{3}(2)));
    end
end

if length(dy(:))==1
    dy{2}=[dy{1},dy{1},2*dy{1}];
    dy{3}=dy{2};
elseif length(dy(:))==2
    if length(dy{2})==1
        dy{2}(2)=dy{2}(1);
        dy{2}(3)=dy{2}(1);
    elseif length(dy{2})==2
        dy{2}(3) = ceil(0.5*(dy{2}(1) + dy{2}(2)));
    end
    dy{3}=dy{2};
else
    if length(dy{2})==1
        dy{2}(2)=dy{2}(1);
        dy{2}(3)=dy{2}(1);
    elseif length(dy{2})==2
        dy{2}(3) = ceil(0.5*(dy{2}(1) + dy{2}(2)));
    end
    if length(dy{3})==1
        dy{3}(2)=dy{3}(1);
        dy{3}(3)=dy{3}(1);
    elseif length(dy{3})==2
        dy{3}(3) = ceil(0.5*(dy{3}(1) + dy{3}(2)));
    end
end

if size(d2,1)==1
    if length(d2{1,1})==1
        d2{1,1} = [d2{1,1};d2{1,1};d2{1,1}];
    elseif length(d2{1,1})==2
        d2{1,1} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(1)+d2{1,1}(2)];
    end
    
    d2{2,1} = [d2{1,1}(1);d2{1,1}(1);d2{1,1}(1);d2{1,1}(2);d2{1,1}(3)];
    d2{3,1} = d2{2,1};
        
    if size(d2,2)==1
        d2{1,2} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(2);d2{1,1}(2);d2{1,1}(3)];
        d2{1,3} = d2{1,2};
    elseif size(d2,2)==2
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        d2{1,3} = d2{1,2};
    else
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        if length(d2{1,3})==1
            d2{1,3}=[d2{1,1}(1);d2{1,3};d2{1,3};d2{1,3};ceil((d2{1,1}(1)+2*d2{1,3})/3)];
        elseif length(d2{1,3})==2
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(2);d2{1,3}(2);ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==3
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(3);ceil(0.5*(d2{1,3}(2)+d2{1,3}(3)));ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==4
            d2{1,3}=[d2{1,3}(1:4);ceil((d2{1,3}(1)+2*d2{1,3}(4))/3)];
        end
        
    end
    d2{2,2} = [d2{2,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,2}(2:4,1)))];
    d2{3,2} = [d2{3,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,2}(2:4,1)))];
    d2{2,3} = [d2{2,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,3}(2:4,1)))];
    d2{3,3} = [d2{3,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,3}(2:4,1)))];
    
    
elseif size(d2,1)==2
    if length(d2{1,1})==1
        d2{1,1}=[d2{1,1};d2{1,1};d2{1,1}];
    elseif length(d2{1,1})==2
        d2{1,1} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(1)+d2{1,1}(2)];
    end
    
    if length(d2{2,1})==1
        d2{2,1}=[d2{2,1};d2{2,1};d2{2,1};d2{2,1};ceil((2*d2{2,1}+d2{1,1}(2))/3)];
    elseif length(d2{2,1})==2
        d2{2,1}=[d2{2,1}(1);d2{2,1}(1);d2{2,1}(1);d2{2,1}(2);ceil(mean(d2{2,1}))];
    elseif length(d2{2,1})==3
        d2{2,1}=[d2{2,1}(1);d2{2,1}(2);ceil(0.5*(d2{2,1}(1)+d2{2,1}(2)));d2{2,1}(3);ceil(mean(d2{2,1}))];
    elseif length(d2{2,1})==4
        d2{2,1}=[d2{2,1}(1:4);ceil((d2{2,1}(4)+2*d2{2,1}(3))/3)];
    end
    d2{3,1} = d2{2,1};
    
    if size(d2,2)==1
        d2{1,2} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(2);d2{1,1}(2);d2{1,1}(3)];
        d2{1,3} = d2{1,2};
    elseif size(d2,2)==2
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        d2{1,3} = d2{1,2};
    else
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        if length(d2{1,3})==1
            d2{1,3}=[d2{1,1}(1);d2{1,3};d2{1,3};d2{1,3};ceil((d2{1,1}(1)+2*d2{1,3})/3)];
        elseif length(d2{1,3})==2
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(2);d2{1,3}(2);ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==3
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(3);ceil(0.5*(d2{1,3}(2)+d2{1,3}(3)));ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==4
            d2{1,3}=[d2{1,3}(1:4);ceil((d2{1,3}(1)+2*d2{1,3}(4))/3)];
        end
    end
        
    %d2{2,2} = [d2{2,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,2}(2:4,1)))];
    d2{3,2} = [d2{3,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,2}(2:4,1)))];
    %d2{2,3} = [d2{2,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,3}(2:4,1)))];
    d2{3,3} = [d2{3,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,3}(2:4,1)))];
    
else
    if length(d2{1,1})==1
        d2{1,1}=[d2{1,1};d2{1,1};d2{1,1}];
    elseif length(d2{1,1})==2
        d2{1,1} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(1)+d2{1,1}(2)];
    end
    
    if length(d2{2,1})==1
        d2{2,1}=[d2{2,1};d2{2,1};d2{2,1};d2{2,1};ceil((2*d2{2,1}+d2{1,1}(2))/3)];
    elseif length(d2{2,1})==2
        d2{2,1}=[d2{2,1}(1);d2{2,1}(1);d2{2,1}(1);d2{2,1}(2);ceil(mean(d2{2,1}))];
    elseif length(d2{2,1})==3
        d2{2,1}=[d2{2,1}(1);d2{2,1}(2);ceil(0.5*(d2{2,1}(1)+d2{2,1}(2)));d2{2,1}(3);ceil(mean(d2{2,1}))];
    elseif length(d2{2,1})==4
        d2{2,1}=[d2{2,1}(1:4);ceil((d2{2,1}(4)+2*d2{2,1}(3))/3)];
    end
    if length(d2{3,1})==1
        d2{3,1}=[d2{3,1};d2{3,1};d2{3,1};d2{3,1};ceil((2*d2{3,1}+d2{1,1}(2))/3)];
    elseif length(d2{3,1})==2
        d2{3,1}=[d2{3,1}(1);d2{3,1}(1);d2{3,1}(1);d2{3,1}(2);ceil(mean(d2{3,1}))];
    elseif length(d2{3,1})==3
        d2{3,1}=[d2{3,1}(1);d2{3,1}(2);ceil(0.5*(d2{3,1}(1)+d2{3,1}(2)));d2{3,1}(3);ceil(mean(d2{3,1}))];
    elseif length(d2{3,1})==4
        d2{3,1}=[d2{3,1}(1:4);ceil((d2{3,1}(4)+2*d2{3,1}(3))/3)];
    end
    
    if size(d2,2)==1
        d2{1,2} = [d2{1,1}(1);d2{1,1}(2);d2{1,1}(2);d2{1,1}(2);d2{1,1}(3)];
        d2{1,3} = d2{1,2};
    elseif size(d2,2)==2
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        d2{1,3} = d2{1,2};
    else
        if length(d2{1,2})==1
            d2{1,2}=[d2{1,1}(1);d2{1,2};d2{1,2};d2{1,2};ceil((d2{1,1}(1)+2*d2{1,2})/3)];
        elseif length(d2{1,2})==2
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(2);d2{1,2}(2);ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==3
            d2{1,2}=[d2{1,2}(1);d2{1,2}(2);d2{1,2}(3);ceil(0.5*(d2{1,2}(2)+d2{1,2}(3)));ceil(mean(d2{1,2}))];
        elseif length(d2{1,2})==4
            d2{1,2}=[d2{1,2}(1:4);ceil((d2{1,2}(1)+2*d2{1,2}(4))/3)];
        end
        if length(d2{1,3})==1
            d2{1,3}=[d2{1,1}(1);d2{1,3};d2{1,3};d2{1,3};ceil((d2{1,1}(1)+2*d2{1,3})/3)];
        elseif length(d2{1,3})==2
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(2);d2{1,3}(2);ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==3
            d2{1,3}=[d2{1,3}(1);d2{1,3}(2);d2{1,3}(3);ceil(0.5*(d2{1,3}(2)+d2{1,3}(3)));ceil(mean(d2{1,3}))];
        elseif length(d2{1,3})==4
            d2{1,3}=[d2{1,3}(1:4);ceil((d2{1,3}(1)+2*d2{1,3}(4))/3)];
        end
    end    
    %d2{2,2} = [d2{2,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,2}(2:4,1)))];
    %d2{3,2} = [d2{3,1}(1:3,1),d2{1,2}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,2}(2:4,1)))];
    %d2{2,3} = [d2{2,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{2,1}(1:3,1)+d2{1,3}(2:4,1)))];
    %d2{3,3} = [d2{3,1}(1:3,1),d2{1,3}(2:4,1),ceil(0.5*(d2{3,1}(1:3,1)+d2{1,3}(2:4,1)))];
    
end

% Specify the spatial domain (s,th \in [I(1),I(2)]
if any(size(I)~=[2,2])
    error('I must be a 2x2 array')
end

% Extract the size of the object: P\in R^(n0 x n0), Rxx\in L_2[x]^(nx x nx),
% Ryy\in L_2[y]^(ny x ny), R22\in L_2[x,y]^(n2 x n2),
n0=n(1);    nx=n(2);    ny=n(3);    n2=n(4);
if n0==0 && nx==0 && ny==0 && n2==0
    error('Error in poslpivar: All dimensions are zero')
end


% % To reduce complexity, allow certain terms to be excluded

% Sorting
% excludeL is a length-16 binary vector of terms to exclude
excludeL = options.exclude;

% In separable case 1, set Ti3 = Ti4 and Zxa=Zxb (so that Rxx{2}=Rxx{3})
if options.sep(1)==1 
    excludeL(4) = 1;
end
% In separable case 2, set Ti6 = Ti7 and Zya=Zyb (so that Ryy{2}=Ryy{3})
if options.sep(2)==1 
    excludeL(7) = 1;
end
% In separable case 3, set Ti9 = Ti10 and Z2ao=Z2bo (so that R22{2,1}=R22{3,1})
if options.sep(3)==1 
    excludeL(10) = 1;
end
% In separable case 4, set Ti11 = Ti12 and Z2oa=Z2ob (so that R22{1,2}=R22{1,3})
if options.sep(4)==1 
    excludeL(12) = 1;
end
% In separable case 5, set Ti13 = Ti14 and Z2aa=Z2ba, and Ti15 = Ti14 and Z2ab=Z2bb  (so that R22{2,2}=R22{3,2} and R22{2,3}=R22{3,3})
% In separable case 6, set Ti13 = Ti15 and Z2aa=Z2ab, and Ti14 = Ti16 and Z2ba=Z2bb  (so that R22{2,2}=R22{2,3} and R22{3,2}=R22{3,3})
if options.sep(5)==1 && options.sep(6) 
    excludeL(14:16) = 1;
elseif options.sep(5)==1
    excludeL([14,16]) = 1;
elseif options.sep(6)==1
    excludeL([15,16]) = 1;
end

if n0==0
    excludeL(1)=1;
end
if nx==0
    excludeL(2:4)=[1 1 1];
end
if ny==0
    excludeL(5:7) = [1 1 1];
end
if n2==0
    excludeL(8:16) = [1 1 1 1 1 1 1 1 1];
end
    
if all(excludeL)
    error('You''re creating an empty dopvar! Please change options.exclude and/or options.sep and try again.')
end


% % Define the variables and multiplier function

% Define the primary variables
pvar ss1 ss2 tt1 tt2;
var11 = ss1; var12 = tt1;
var21 = ss2; var22 = tt2;
var1 = [var11;var21];
var2 = [var12;var22];

% Define dummy variables for integration (replace eta)
rr1=polynomial(1,1,{'rr1'},[1 1]);
rr2=polynomial(1,1,{'rr2'},[1 1]);
rr = [rr1;rr2];

% Define the multiplier function to be used later

psatz = abs(options.psatz);
if psatz==0
    gss=polynomial(1);
    
elseif psatz==1
    gss=(var11-I(1,1))*(I(1,2)-var11)*(var21-I(2,1))*(I(2,2)-var21);
    
elseif psatz==2
    cntr = [mean(I(1,:)),mean(I(2,:))];
    rds = norm([I(1,2),I(2,2)] - cntr);
    
    gss = (rds^2 - (var11-cntr(1))^2 - (var21-cntr(2))^2);
    
else
    error('options.psatz can only assume values 0, 1, and 2')
end

% ONLY FOR TESTING PURPOSES: Take square of psatz polynomial
% if sign(options.psatz)==-1
%     gss = gss^2;  
% end

% Introduce varitations with adjusted variables for later purposes
gtt = subs(gss,var1,var2);
grr = subs(gss,var1,rr);
gts = subs(gss,var11,var12);
gst = subs(gss,var21,var22);
grs = subs(gss,var11,rr1);
gsr = subs(gss,var21,rr2);
grt = subs(gtt,var12,rr1);
gtr = subs(gtt,var22,rr2);

% Introduce integrated versions for later purposes
gis = int(gss,var11,I(1,1),I(1,2));
gsi = int(gss,var21,I(2,1),I(2,2));
git = subs(gis,var21,var22);
gti = subs(gsi,var11,var12);
gir = subs(gis,var21,rr2);
gri = subs(gsi,var11,rr1);

gii = int(gis,var21,I(2,1),I(2,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Construct the monomial bases Z % % %

% % Build the monomials

    % Constructing Zxo(ss1)
nZxo = dx{1}+1;
Zxo_degmat = (0:dx{1})';
Zxo_coeff = speye(nZxo);
Zxo_varname = var11.varname;
Zxo_matdim = [nZxo 1];
Zxo_s = polynomial(Zxo_coeff,Zxo_degmat,Zxo_varname,Zxo_matdim);


    % Constructing Zxa(ss1,tt1) and Zxb(ss1,tt1)
% In this implementation, Zxa will have degree dx{2}(2) in tt1 and degree
% d{2}(1) in ss1 and max degree of ss1+tt1 is d{2}(3). Similarly for Zxb(ss1,tt1)

Zxa_degmat = [repmat((0:dx{2}(1))',dx{2}(2)+1,1),vec(repmat(0:dx{2}(2),dx{2}(1)+1,1))];
Zxa_degmat(sum(Zxa_degmat,2)>dx{2}(3),:)= [];
nZxa=size(Zxa_degmat,1);
Zxa_coeff = speye(nZxa);
Zxa_varname = [var11.varname; var12.varname];
Zxa_matdim = [nZxa 1];
Zxa_st=polynomial(Zxa_coeff,Zxa_degmat,Zxa_varname,Zxa_matdim);

Zxb_degmat = [repmat((0:dx{3}(1))',dx{3}(2)+1,1),vec(repmat(0:dx{3}(2),dx{3}(1)+1,1))];
Zxb_degmat(sum(Zxb_degmat,2)>dx{3}(3),:)= [];
nZxb=size(Zxb_degmat,1);
Zxb_coeff = speye(nZxb);
Zxb_varname = [var11.varname; var12.varname];
Zxb_matdim = [nZxb 1];
Zxb_st = polynomial(Zxb_coeff,Zxb_degmat,Zxb_varname,Zxb_matdim);


% % % % %

    % Constructing Zyo(ss2)
nZyo = dy{1}+1;
Zyo_degmat = (0:dy{1})';
Zyo_coeff = speye(nZyo);
Zyo_varname = var21.varname;
Zyo_matdim = [nZyo 1];
Zyo_s = polynomial(Zyo_coeff,Zyo_degmat,Zyo_varname,Zyo_matdim);

    % Constructing Zya(ss2,tt2) and Zyb(ss2,tt2)
% In this implementation, Zxa will have degree dx{2}(2) in tt1 and degree
% d{2}(1) in ss1 and max degree of ss1+tt1 is d{2}(3). Similarly for Zxb(ss1,tt1)

Zya_degmat = [repmat((0:dy{2}(1))',dy{2}(2)+1,1),vec(repmat(0:dy{2}(2),dy{2}(1)+1,1))];
Zya_degmat(sum(Zya_degmat,2)>dy{2}(3),:) = [];
nZya = size(Zya_degmat,1);
Zya_coeff = speye(nZya);
Zya_varname = [var21.varname; var22.varname];
Zya_matdim = [nZya 1];
Zya_st = polynomial(Zya_coeff,Zya_degmat,Zya_varname,Zya_matdim);

Zyb_degmat = [repmat((0:dy{3}(1))',dy{3}(2)+1,1),vec(repmat(0:dy{3}(2),dy{3}(1)+1,1))];
Zyb_degmat(sum(Zyb_degmat,2)>dy{3}(3),:) = [];
nZyb = size(Zyb_degmat,1);
Zyb_coeff = speye(nZyb);
Zyb_varname = [var21.varname; var22.varname];
Zyb_matdim = [nZyb 1];
Zyb_st = polynomial(Zyb_coeff,Zyb_degmat,Zyb_varname,Zyb_matdim);


% % % % %

    % Constructing Z2oo(ss1,ss2)
Z2oo_degmat = [repmat((0:d2{1}(1))',d2{1}(2)+1,1),vec(repmat(0:d2{1}(2),d2{1}(1)+1,1))];
Z2oo_degmat(sum(Z2oo_degmat,2)>d2{1}(3),:) = [];
nZ2oo = size(Z2oo_degmat,1);
Z2oo_coeff = speye(nZ2oo);
Z2oo_varname = [var11.varname; var21.varname];
Z2oo_matdim = [nZ2oo 1];
Z2oo_ss = polynomial(Z2oo_coeff,Z2oo_degmat,Z2oo_varname,Z2oo_matdim);

% %
    % Constructing Z2ao(ss1,ss2,tt1) and Z2bo(ss1,ss2,tt1)
Z2ao_degmat = [repmat((0:d2{2,1}(1))',(d2{2,1}(2)+1)*(d2{2,1}(4)+1),1),...
                repmat(vec(repmat(0:d2{2,1}(2),d2{2,1}(1)+1,1)),(d2{2,1}(4)+1),1),...
                vec(repmat(0:d2{2,1}(4),(d2{2,1}(1)+1)*(d2{2,1}(2)+1),1))];
Z2ao_degmat(sum(Z2ao_degmat(:,[1,2]),2)>d2{2,1}(3),:) = [];
Z2ao_degmat(sum(Z2ao_degmat,2)>d2{2,1}(5),:) = [];
nZ2ao = size(Z2ao_degmat,1);
Z2ao_coeff = speye(nZ2ao);
Z2ao_varname = [var11.varname; var12.varname; var21.varname];
Z2ao_matdim = [nZ2ao 1];
Z2ao_sst = polynomial(Z2ao_coeff,Z2ao_degmat,Z2ao_varname,Z2ao_matdim);

Z2bo_degmat = [repmat((0:d2{3,1}(1))',(d2{3,1}(2)+1)*(d2{3,1}(4)+1),1),...
                repmat(vec(repmat(0:d2{3,1}(2),d2{3,1}(1)+1,1)),(d2{3,1}(4)+1),1),...
                vec(repmat(0:d2{3,1}(4),(d2{3,1}(1)+1)*(d2{3,1}(2)+1),1))];
Z2bo_degmat(sum(Z2bo_degmat(:,[1,2]),2)>d2{3,1}(3),:) = [];
Z2bo_degmat(sum(Z2bo_degmat,2)>d2{3,1}(5),:) = [];
nZ2bo = size(Z2bo_degmat,1);
Z2bo_coeff = speye(nZ2bo);
Z2bo_varname = [var11.varname; var12.varname; var21.varname];
Z2bo_matdim = [nZ2bo 1];
Z2bo_sst = polynomial(Z2bo_coeff,Z2bo_degmat,Z2bo_varname,Z2bo_matdim);

%
    % Constructing Z2oa(ss1,ss2,tt2) and Z2ob(ss1,ss2,tt2)
Z2oa_degmat = [repmat((0:d2{1,2}(1))',(d2{1,2}(2)+1)*(d2{1,2}(3)+1),1),...
                repmat(vec(repmat(0:d2{1,2}(2),d2{1,2}(1)+1,1)),(d2{1,2}(3)+1),1),...
                vec(repmat(0:d2{1,2}(3),(d2{1,2}(1)+1)*(d2{1,2}(2)+1),1))];
Z2oa_degmat(sum(Z2oa_degmat(:,[2,3]),2)>d2{1,2}(4),:) = [];
Z2oa_degmat(sum(Z2oa_degmat,2)>d2{1,2}(5),:) = [];
nZ2oa = size(Z2oa_degmat,1);
Z2oa_coeff = speye(nZ2oa);
Z2oa_varname = [var11.varname; var21.varname; var22.varname];
Z2oa_matdim = [nZ2oa 1];
Z2oa_sst = polynomial(Z2oa_coeff,Z2oa_degmat,Z2oa_varname,Z2oa_matdim);

Z2ob_degmat = [repmat((0:d2{1,3}(1))',(d2{1,3}(2)+1)*(d2{1,3}(3)+1),1),...
                repmat(vec(repmat(0:d2{1,3}(2),d2{1,3}(1)+1,1)),(d2{1,3}(3)+1),1),...
                vec(repmat(0:d2{1,3}(3),(d2{1,3}(1)+1)*(d2{1,3}(2)+1),1))];
Z2ob_degmat(sum(Z2ob_degmat(:,[2,3]),2)>d2{1,3}(4),:) = [];
Z2ob_degmat(sum(Z2ob_degmat,2)>d2{1,3}(5),:) = [];
nZ2ob = size(Z2ob_degmat,1);
Z2ob_coeff = speye(nZ2ob);
Z2ob_varname = [var11.varname; var21.varname; var22.varname];
Z2ob_matdim = [nZ2ob 1];
Z2ob_sst = polynomial(Z2ob_coeff,Z2ob_degmat,Z2ob_varname,Z2ob_matdim);

%
    % Constructing Z2aa(ss1,ss2,tt1,tt2) ... Z2bb(ss1,ss2,tt1,tt2)
Z2aa_degmat = [repmat((0:d2{2,2}(1,1))',(d2{2,2}(2,1)+1)*(d2{2,2}(1,2)+1)*(d2{2,2}(2,2)+1),1),...
                repmat(vec(repmat(0:d2{2,2}(2,1),d2{2,2}(1,1)+1,1)),(d2{2,2}(1,2)+1)*(d2{2,2}(2,2)+1),1),...
                repmat(vec(repmat(0:d2{2,2}(1,2),(d2{2,2}(1,1)+1)*(d2{2,2}(2,1)+1),1)),(d2{2,2}(2,2)+1),1),...
                vec(repmat(0:d2{2,2}(2,2),(d2{2,2}(1,1)+1)*(d2{2,2}(2,1)+1)*(d2{2,2}(1,2)+1),1))];
Z2aa_degmat(sum(Z2aa_degmat(:,[1,2]),2)>d2{2,2}(3,1),:) = [];
Z2aa_degmat(sum(Z2aa_degmat(:,[3,4]),2)>d2{2,2}(3,2),:) = [];
Z2aa_degmat(sum(Z2aa_degmat(:,[1,3]),2)>d2{2,2}(1,3),:) = [];
Z2aa_degmat(sum(Z2aa_degmat(:,[2,4]),2)>d2{2,2}(2,3),:) = [];
Z2aa_degmat(sum(Z2aa_degmat,2)>d2{2,2}(3,3),:) = [];
nZ2aa = size(Z2aa_degmat,1);
Z2aa_coeff = speye(nZ2aa);
Z2aa_varname = [var11.varname; var12.varname; var21.varname; var22.varname];
Z2aa_matdim = [nZ2aa 1];
Z2aa_sstt = polynomial(Z2aa_coeff,Z2aa_degmat,Z2aa_varname,Z2aa_matdim);

Z2ba_degmat = [repmat((0:d2{3,2}(1,1))',(d2{3,2}(2,1)+1)*(d2{3,2}(1,2)+1)*(d2{3,2}(2,2)+1),1),...
                repmat(vec(repmat(0:d2{3,2}(2,1),d2{3,2}(1,1)+1,1)),(d2{3,2}(1,2)+1)*(d2{3,2}(2,2)+1),1),...
                repmat(vec(repmat(0:d2{3,2}(1,2),(d2{3,2}(1,1)+1)*(d2{3,2}(2,1)+1),1)),(d2{3,2}(2,2)+1),1),...
                vec(repmat(0:d2{3,2}(2,2),(d2{3,2}(1,1)+1)*(d2{3,2}(2,1)+1)*(d2{3,2}(1,2)+1),1))];
Z2ba_degmat(sum(Z2ba_degmat(:,[1,2]),2)>d2{3,2}(3,1),:) = [];
Z2ba_degmat(sum(Z2ba_degmat(:,[3,4]),2)>d2{3,2}(3,2),:) = [];
Z2ba_degmat(sum(Z2ba_degmat(:,[1,3]),2)>d2{3,2}(1,3),:) = [];
Z2ba_degmat(sum(Z2ba_degmat(:,[2,4]),2)>d2{3,2}(2,3),:) = [];
Z2ba_degmat(sum(Z2ba_degmat,2)>d2{3,2}(3,3),:) = [];
nZ2ba = size(Z2ba_degmat,1);
Z2ba_coeff = speye(nZ2ba);
Z2ba_varname = [var11.varname; var12.varname; var21.varname; var22.varname];
Z2ba_matdim = [nZ2ba 1];
Z2ba_sstt = polynomial(Z2ba_coeff,Z2ba_degmat,Z2ba_varname,Z2ba_matdim);

Z2ab_degmat = [repmat((0:d2{2,3}(1,1))',(d2{2,3}(2,1)+1)*(d2{2,3}(1,2)+1)*(d2{2,3}(2,2)+1),1),...
                repmat(vec(repmat(0:d2{2,3}(2,1),d2{2,3}(1,1)+1,1)),(d2{2,3}(1,2)+1)*(d2{2,3}(2,2)+1),1),...
                repmat(vec(repmat(0:d2{2,3}(1,2),(d2{2,3}(1,1)+1)*(d2{2,3}(2,1)+1),1)),(d2{2,3}(2,2)+1),1),...
                vec(repmat(0:d2{2,3}(2,2),(d2{2,3}(1,1)+1)*(d2{2,3}(2,1)+1)*(d2{2,3}(1,2)+1),1))];
Z2ab_degmat(sum(Z2ab_degmat(:,[1,2]),2)>d2{2,3}(3,1),:) = [];
Z2ab_degmat(sum(Z2ab_degmat(:,[3,4]),2)>d2{2,3}(3,2),:) = [];
Z2ab_degmat(sum(Z2ab_degmat(:,[1,3]),2)>d2{2,3}(1,3),:) = [];
Z2ab_degmat(sum(Z2ab_degmat(:,[2,4]),2)>d2{2,3}(2,3),:) = [];
Z2ab_degmat(sum(Z2ab_degmat,2)>d2{2,3}(3,3),:) = [];
nZ2ab = size(Z2ab_degmat,1);
Z2ab_coeff = speye(nZ2ab);
Z2ab_varname = [var11.varname; var12.varname; var21.varname; var22.varname];
Z2ab_matdim = [nZ2ab 1];
Z2ab_sstt = polynomial(Z2ab_coeff,Z2ab_degmat,Z2ab_varname,Z2ab_matdim);

Z2bb_degmat = [repmat((0:d2{3,3}(1,1))',(d2{3,3}(2,1)+1)*(d2{3,3}(1,2)+1)*(d2{3,3}(2,2)+1),1),...
                repmat(vec(repmat(0:d2{3,3}(2,1),d2{3,3}(1,1)+1,1)),(d2{3,3}(1,2)+1)*(d2{3,3}(2,2)+1),1),...
                repmat(vec(repmat(0:d2{3,3}(1,2),(d2{3,3}(1,1)+1)*(d2{3,3}(2,1)+1),1)),(d2{3,3}(2,2)+1),1),...
                vec(repmat(0:d2{3,3}(2,2),(d2{3,3}(1,1)+1)*(d2{3,3}(2,1)+1)*(d2{3,3}(1,2)+1),1))];
Z2bb_degmat(sum(Z2bb_degmat(:,[1,2]),2)>d2{3,3}(3,1),:) = [];
Z2bb_degmat(sum(Z2bb_degmat(:,[3,4]),2)>d2{3,3}(3,2),:) = [];
Z2bb_degmat(sum(Z2bb_degmat(:,[1,3]),2)>d2{3,3}(1,3),:) = [];
Z2bb_degmat(sum(Z2bb_degmat(:,[2,4]),2)>d2{3,3}(2,3),:) = [];
Z2bb_degmat(sum(Z2bb_degmat,2)>d2{3,3}(3,3),:) = [];
nZ2bb = size(Z2bb_degmat,1);
Z2bb_coeff = speye(nZ2bb);
Z2bb_varname = [var11.varname; var12.varname; var21.varname; var22.varname];
Z2bb_matdim = [nZ2bb 1];
Z2bb_sstt = polynomial(Z2bb_coeff,Z2bb_degmat,Z2bb_varname,Z2bb_matdim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Construct the left and right monomial bases ZL and ZR % % %

% Introduce monomials with adjusted variables
Zxo_t = subs(Zxo_s,var11,var12);
Zxa_rt = subs(Zxa_st,var11,rr1);
Zxa_rs = subs(Zxa_rt,var12,var11);
Zxb_rt = subs(Zxb_st,var11,rr1);
Zxb_rs = subs(Zxb_rt,var12,var11);

Zyo_t = subs(Zyo_s,var21,var22);
Zya_rt = subs(Zya_st,var21,rr2);
Zya_rs = subs(Zya_rt,var22,var21);
Zyb_rt = subs(Zyb_st,var21,rr2);
Zyb_rs = subs(Zyb_rt,var22,var21);

Z2oo_tt = subs(Z2oo_ss,var1,var2);

Z2ao_rst = subs(Z2ao_sst,var11,rr1);
Z2ao_rtt = subs(Z2ao_rst,var21,var22);
Z2ao_rss = subs(Z2ao_rst,var12,var11);
Z2bo_rst = subs(Z2bo_sst,var11,rr1);
Z2bo_rtt = subs(Z2bo_rst,var21,var22);
Z2bo_rss = subs(Z2bo_rst,var12,var11);

Z2oa_srt = subs(Z2oa_sst,var21,rr2);
Z2oa_trt = subs(Z2oa_srt,var11,var12);
Z2oa_srs = subs(Z2oa_srt,var22,var21);
Z2ob_srt = subs(Z2ob_sst,var21,rr2);
Z2ob_trt = subs(Z2ob_srt,var11,var12);
Z2ob_srs = subs(Z2ob_srt,var22,var21);

Z2aa_rrtt = subs(Z2aa_sstt,var1,rr);
Z2aa_rrss = subs(Z2aa_rrtt,var2,var1);
Z2ba_rrtt = subs(Z2ba_sstt,var1,rr);
Z2ba_rrss = subs(Z2ba_rrtt,var2,var1);
Z2ab_rrtt = subs(Z2ab_sstt,var1,rr);
Z2ab_rrss = subs(Z2ab_rrtt,var2,var1);
Z2bb_rrtt = subs(Z2bb_sstt,var1,rr);
Z2bb_rrss = subs(Z2bb_rrtt,var2,var1);


% Collect the different monomials in cells ZL and ZR
includeL = ~excludeL;
ZL = cell(1,sum(includeL));
ZR = cell(1,sum(includeL));

mdim = [];
ndim = [];
indx = 1;
if includeL(1)
    ZL{indx} = 1;
    ZR{indx} = 1;
    mdim = [mdim;n0];
    ndim = [ndim;n0];
    indx = indx+1;
end
if includeL(2)
    ZL{indx} = Zxo_s;
    ZR{indx} = Zxo_t;
    mdim = [mdim;nx];
    ndim = [ndim;nx];
    indx = indx+1;
end
if includeL(3)
    ZL{indx} = Zxa_rs;
    ZR{indx} = Zxa_rt;
    mdim = [mdim;nx];
    ndim = [ndim;nx];
    indx = indx+1;
end
if includeL(4)
    ZL{indx} = Zxb_rs;
    ZR{indx} = Zxb_rt;
    mdim = [mdim;nx];
    ndim = [ndim;nx];
    indx = indx+1;
end
if includeL(5)
    ZL{indx} = Zyo_s;
    ZR{indx} = Zyo_t;
    mdim = [mdim;ny];
    ndim = [ndim;ny];
    indx = indx+1;
end
if includeL(6)
    ZL{indx} = Zya_rs;
    ZR{indx} = Zya_rt;
    mdim = [mdim;ny];
    ndim = [ndim;ny];
    indx = indx+1;
end
if includeL(7)
    ZL{indx} = Zyb_rs;
    ZR{indx} = Zyb_rt;
    mdim = [mdim;ny];
    ndim = [ndim;ny];
    indx = indx+1;
end
if includeL(8)
    ZL{indx} = Z2oo_ss;
    ZR{indx} = Z2oo_tt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(9)
    ZL{indx} = Z2ao_rss;
    ZR{indx} = Z2ao_rtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(10)
    ZL{indx} = Z2bo_rss;
    ZR{indx} = Z2bo_rtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(11)
    ZL{indx} = Z2oa_srs;
    ZR{indx} = Z2oa_trt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(12)
    ZL{indx} = Z2ob_srs;
    ZR{indx} = Z2ob_trt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(13)
    ZL{indx} = Z2aa_rrss;
    ZR{indx} = Z2aa_rrtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(14)
    ZL{indx} = Z2ba_rrss;
    ZR{indx} = Z2ba_rrtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(15)
    ZL{indx} = Z2ab_rrss;
    ZR{indx} = Z2ab_rrtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    indx = indx+1;
end
if includeL(16)
    ZL{indx} = Z2bb_rrss;
    ZR{indx} = Z2bb_rrtt;
    mdim = [mdim;n2];
    ndim = [ndim;n2];
    %indx = indx+1;
end

% Compute the product N = ZL'*T*ZR for positive decision matrix T
[prog,N] = sosquadvar(prog,ZL,ZR,mdim,ndim,'pos');
ij = cumsum(includeL);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Build the positive decision operator % % %

dopvar2d Pop;
Pop.dim = [n0,n0; nx,nx; ny,ny; n2,n2];
Pop.var1 = [var11;var21];
Pop.var2 = [var12;var22];
Pop.I = I;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Build R00, R0x, R0y and R02
if includeL(1)
    %Pop.R00 = Pop.R00 + gii*Q{1,1};
    Pop.R00 = Pop.R00 + gii*N{ij(1),ij(1)};
    if includeL(2)
        %Pop.R0x = Pop.R0x + gsi*Q{1,2}*bZxo_s;
        Pop.R0x = Pop.R0x + gsi*subs(N{ij(1),ij(2)},var12,var11);
    end
    if includeL(3) && ~options.sep(1)==1
        %Pop.R0x = Pop.R0x + int(gti*Q{1,3}*bZxa_ts,var12,var11,I(1,2));
        Pop.R0x = Pop.R0x + int(gri*subs(N{ij(1),ij(3)},var12,var11),rr1,var11,I(1,2));
    elseif includeL(3) && options.sep(1)==1
        %Pop.R0x = Pop.R0x + int(gti*Q{1,3}*bZxa_ts,var12,I(1,1),I(1,2));
        Pop.R0x = Pop.R0x + int(gri*subs(N{ij(1),ij(3)},var12,var11),rr1,I(1,1),I(1,2));
    end
    if includeL(4)
        %Pop.R0x = Pop.R0x + int(gti*Q{1,4}*bZxb_ts,var12,I(1,1),var11);
        Pop.R0x = Pop.R0x + int(gri*subs(N{ij(1),ij(4)},var12,var11),rr1,I(1,1),var11);
    end
    
    if includeL(5)
        %Pop.R0y = Pop.R0y + gis*Q{1,5}*bZyo_s;
        Pop.R0y = Pop.R0y + gis*subs(N{ij(1),ij(5)},var22,var21);
    end
    if includeL(6) && ~options.sep(2)==1
        %Pop.R0y = Pop.R0y + int(git*Q{1,6}*bZya_ts,var22,var21,I(2,2));
        Pop.R0y = Pop.R0y + int(gir*subs(N{ij(1),ij(6)},var22,var21),rr2,var21,I(2,2));
    elseif includeL(3) && options.sep(2)==1
        %Pop.R0y = Pop.R0y + int(git*Q{1,6}*bZya_ts,var22,I(2,1),I(2,2));
        Pop.R0y = Pop.R0y + int(gir*subs(N{ij(1),ij(6)},var22,var21),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        %Pop.R0y = Pop.R0y + int(git*Q{1,7}*bZyb_ts,var22,I(2,1),var21);
        Pop.R0y = Pop.R0y + int(gir*subs(N{ij(1),ij(7)},var22,var21),rr2,I(2,1),var21);
    end
    
    if includeL(8)
        %Pop.R02 = Pop.R02 + gss*Q{1,8}*bZ2oo_ss; 
        Pop.R02 = Pop.R02 + gss*subs(N{ij(1),ij(8)},var2,var1);
    end
    if includeL(9) && ~options.sep(3)==1
        %Pop.R02 = Pop.R02 + int(gts*Q{1,9}*bZ2ao_tss,var12,var11,I(1,2));
        Pop.R02 = Pop.R02 + int(grs*subs(N{ij(1),ij(9)},var2,var1),rr1,var11,I(1,2));
    elseif includeL(9) && options.sep(3)==1
        %Pop.R02 = Pop.R02 + int(gts*Q{1,9}*bZ2ao_tss,var12,I(1,1),I(1,2));
        Pop.R02 = Pop.R02 + int(grs*subs(N{ij(1),ij(9)},var2,var1),rr1,I(1,1),I(1,2));
    end
    if includeL(10)
        %Pop.R02 = Pop.R02 + int(gts*Q{1,10}*bZ2bo_tss,var12,I(1,1),var11);
        Pop.R02 = Pop.R02 + int(grs*subs(N{ij(1),ij(10)},var2,var1),rr1,I(1,1),var11);
    end
    if includeL(11) && ~options.sep(4)==1
        %Pop.R02 = Pop.R02 + int(gst*Q{1,11}*bZ2oa_sts,var22,var21,I(2,2));
        Pop.R02 = Pop.R02 + int(gsr*subs(N{ij(1),ij(11)},var2,var1),rr2,var21,I(2,2));
    elseif includeL(11) && options.sep(4)==1
        %Pop.R02 = Pop.R02 + int(gst*Q{1,11}*bZ2oa_sts,var22,I(2,1),I(2,2));
        Pop.R02 = Pop.R02 + int(gsr*subs(N{ij(1),ij(11)},var2,var1),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        %Pop.R02 = Pop.R02 + int(gst*Q{1,12}*bZ2ob_sts,var22,I(2,1),var21);
        Pop.R02 = Pop.R02 + int(gsr*subs(N{ij(1),ij(12)},var2,var1),rr2,I(2,1),var21);
    end
    if includeL(13) && ~options.sep(5)==1 && ~options.sep(6)
        %Pop.R02 = Pop.R02 + int(int(gtt*Q{1,13}*bZ2aa_tsts,var12,var11,I(1,2)),var22,var21,I(2,2));
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(13)},var2,var1),rr1,var11,I(1,2)),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(13)},var2,var1),rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(13)},var2,var1),rr1,var11,I(1,2)),rr2,I(2,1),I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.R02 = Pop.R02 + int(int(gtt*Q{1,13}*bZ2aa_tsts,var12,I(1,1),I(1,2)),var22,I(2,1),I(2,2));
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(13)},var2,var1),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(14) && ~options.sep(6)
        %Pop.R02 = Pop.R02 + int(int(gtt*Q{1,14}*bZ2ba_tsts,var12,I(1,1),var11),var22,var21,I(2,2));
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(14)},var2,var1),rr1,I(1,1),var11),rr2,var21,I(2,2));
    elseif includeL(14) && options.sep(6)
        %
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(14)},var2,var1),rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~options.sep(5)
        %Pop.R02 = Pop.R02 + int(int(gtt*Q{1,15}*bZ2ab_tsts,var12,var11,I(1,2)),var22,I(2,1),var21);
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(15)},var2,var1),rr1,var11,I(1,2)),rr2,I(2,1),var21);
    elseif includeL(15) && options.sep(5)
        %
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(15)},var2,var1),rr1,I(1,1),I(1,2)),rr2,I(2,1),var21);
    end
    if includeL(16) 
        %Pop.R02 = Pop.R02 + int(int(gtt*Q{1,16}*bZ2bb_tsts,var12,I(1,1),var11),var22,I(2,1),var21);
        Pop.R02 = Pop.R02 + int(int(grr*subs(N{ij(1),ij(16)},var2,var1),rr1,I(1,1),var11),rr2,I(2,1),var21);
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Build Rxx, Rxy and Rx2
if includeL(2)
    %Pop.Rxx{1} = constructR0(Q{2,2},Zxo_s,gsi,nx,nZxo); 
    Pop.Rxx{1} = gsi * subs(N{ij(2),ij(2)},var12,var11); 
    if includeL(3) && ~options.sep(1)==1
        %Pop.Rxx{2} = Pop.Rxx{2} + constructA1(Q{2,3},Zxo_s,Zxa_st,gsi,nx,nZxo,nZxa,var11,var12); 
        Pop.Rxx{2} = Pop.Rxx{2} + gsi * subs(N{ij(2),ij(3)},rr1,var11); 
    elseif includeL(3) && options.sep(1)==1
        %Pop.Rxx{2} = Pop.Rxx{2} + constructA1(Q{2,3},Zxo_s,Zxa_st,gsi,nx,nZxo,nZxa,var11,var12) + constructB1(Q{3,2},Zxo_s,Zxa_st,gti,nx,nZxo,nZxa,var11,var12);
        Pop.Rxx{2} = Pop.Rxx{2} + gsi * subs(N{ij(2),ij(3)},rr1,var11) + gti * subs(N{ij(3),ij(2)},rr1,var12); 
    end
    if includeL(4)
        %Pop.Rxx{2} = Pop.Rxx{2} + constructB1(Q{4,2},Zxo_s,Zxa_st,gti,nx,nZxo,nZxa,var11,var12); 
        Pop.Rxx{2} = Pop.Rxx{2} + gti * subs(N{ij(4),ij(2)},rr1,var12); 
    end
    
    if includeL(5)
        %Pop.Rxy = Pop.Rxy + constructPxy(Q{2,5},Zxo_s,Zyo_s,gss,nx,ny,nZxo,nZyo);
        Pop.Rxy = Pop.Rxy + gss * subs(N{ij(2),ij(5)},var22,var21);
    end
    if includeL(6) && ~options.sep(2)==1
        %Pop.Rxy = Pop.Rxy + int(constructPxy(Q{2,6},Zxo_s,Zya_st,gst,nx,ny,nZxo,nZya),var22,var21,I(2,2));
        Pop.Rxy = Pop.Rxy + int(gsr * subs(N{ij(2),ij(6)},var22,var21),rr2,var21,I(2,2));
    elseif includeL(6) && options.sep(2)==1
        %Pop.Rxy = Pop.Rxy + int(constructPxy(Q{2,6},Zxo_s,Zya_st,gst,nx,ny,nZxo,nZya),var22,I(2,1),I(2,2));
        Pop.Rxy = Pop.Rxy + int(gsr * subs(N{ij(2),ij(6)},var22,var21),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        %Pop.Rxy = Pop.Rxy + int(constructPxy(Q{2,7},Zxo_s,Zyb_st,gst,nx,ny,nZxo,nZyb),var22,I(2,1),var21);
        Pop.Rxy = Pop.Rxy + int(gsr * subs(N{ij(2),ij(7)},var22,var21),rr2,I(2,1),var21);
    end
    
    if includeL(8)
        %Pop.Rx2{1} = Pop.Rx2{1} + constructPx2a(Q{2,8},Zxo_s,Z2oo_ss,gss,nx,n2,nZxo,nZ2oo,[0;0]);
        Pop.Rx2{1} = Pop.Rx2{1} + gss * subs(N{ij(2),ij(8)},var2,var1);
    end
    if includeL(11) && ~options.sep(4)==1
        %Pop.Rx2{1} = Pop.Rx2{1} + int(constructPx2a(Q{2,11},Zxo_s,Z2oa_sst,gst,nx,n2,nZxo,nZ2oa,[0;1]),var22,var21,I(2,2));
        Pop.Rx2{1} = Pop.Rx2{1} + int(gsr * subs(N{ij(2),ij(11)},var2,var1),rr2,var21,I(2,2));
    elseif includeL(11) && options.sep(4)==1
        %Pop.Rx2{1} = Pop.Rx2{1} + int(constructPx2a(Q{2,11},Zxo_s,Z2oa_sst,gst,nx,n2,nZxo,nZ2oa,[0;1]),var22,I(2,1),I(2,2));
        Pop.Rx2{1} = Pop.Rx2{1} + int(gsr * subs(N{ij(2),ij(11)},var2,var1),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        %Pop.Rx2{1} = Pop.Rx2{1} + int(constructPx2a(Q{2,12},Zxo_s,Z2ob_sst,gst,nx,n2,nZxo,nZ2ob,[0;1]),var22,I(2,1),var21);
        Pop.Rx2{1} = Pop.Rx2{1} + int(gsr * subs(N{ij(2),ij(12)},var2,var1),rr2,I(2,1),var21);
    end
    
    if includeL(9) && ~options.sep(3)
        %Pop.Rx2{2} = Pop.Rx2{2} + constructPx2a(Q{2,9},Zxo_s,Z2ao_sst,gss,nx,n2,nZxo,nZ2ao,[0;0]);
        Pop.Rx2{2} = Pop.Rx2{2} + gss * subs(N{ij(2),ij(9)},[rr1;var22],var1);
    elseif includeL(9) && options.sep(3)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + gss * subs(N{ij(2),ij(9)},[rr1;var22],var1);
        Pop.Rx2{3} = Pop.Rx2{3} + gss * subs(N{ij(2),ij(9)},[rr1;var22],var1);
    end
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2a(Q{2,13},Zxo_s,Z2aa_sstt,gst,nx,n2,nZxo,nZ2aa,[0;1]),var22,var21,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(gsr * subs(N{ij(2),ij(13)},[rr1;var22],var1),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(gsr * subs(N{ij(2),ij(13)},[rr1;var22],var1),rr2,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(gsr * subs(N{ij(2),ij(13)},[rr1;var22],var1),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2a(Q{2,13},Zxo_s,Z2aa_sstt,gst,nx,n2,nZxo,nZ2aa,[0;1]),var22,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(gsr * subs(N{ij(2),ij(13)},[rr1;var22],var1),rr2,I(2,1),I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(gsr * subs(N{ij(2),ij(13)},[rr1;var22],var1),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(gsr * subs(N{ij(2),ij(13)},[rr1;var22],var1),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~options.sep(6)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2a(Q{2,15},Zxo_s,Z2ab_sstt,gst,nx,n2,nZxo,nZ2ab,[0;1]),var22,I(2,1),var21);
        Pop.Rx2{2} = Pop.Rx2{2} + int(gsr * subs(N{ij(2),ij(15)},[rr1;var22],var1),rr2,I(2,1),var21);
    elseif includeL(15) && options.sep(6)
        Pop.Rx2{2} = Pop.Rx2{2} + int(gsr * subs(N{ij(2),ij(15)},[rr1;var22],var1),rr2,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(gsr * subs(N{ij(2),ij(15)},[rr1;var22],var1),rr2,I(2,1),var21);
    end
    
    if includeL(10)
        %Pop.Rx2{3} = Pop.Rx2{3} + constructPx2a(Q{2,10},Zxo_s,Z2bo_sst,gss,nx,n2,nZxo,nZ2bo,[0;0]);
        Pop.Rx2{3} = Pop.Rx2{3} + gss * subs(N{ij(2),ij(10)},[rr1;var22],var1);
    end
    if includeL(14) && ~options.sep(5)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2a(Q{2,14},Zxo_s,Z2ba_sstt,gst,nx,n2,nZxo,nZ2ba,[0;1]),var22,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(gsr * subs(N{ij(2),ij(14)},[rr1;var22],var1),rr2,var21,I(2,2));
    elseif includeL(14) && options.sep(5)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2a(Q{2,14},Zxo_s,Z2ba_sstt,gst,nx,n2,nZxo,nZ2ba,[0;1]),var22,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(gsr * subs(N{ij(2),ij(14)},[rr1;var22],var1),rr2,I(2,1),I(2,2));
    end
    if includeL(16)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2a(Q{2,16},Zxo_s,Z2bb_sstt,gst,nx,n2,nZxo,nZ2bb,[0;1]),var22,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(gsr * subs(N{ij(2),ij(16)},[rr1;var22],var1),rr2,I(2,1),var21);
    end
end

if includeL(3) && ~options.sep(1)==1
    %Pop.Rxx{2} = Pop.Rxx{2} + int(constructCA1(Q{3,3},Zxa_st,gri,nx,nZxa,var11,var12,rr1),rr1,var11,I(1,2));
    Pop.Rxx{2} = Pop.Rxx{2} + int(gri * N{ij(3),ij(3)},rr1,var11,I(1,2));
    if includeL(4)
        %Pop.Rxx{2} = Pop.Rxx{2} + int(constructCA1(Q{4,3},Zxa_st,gri,nx,nZxa,var11,var12,rr1),rr1,var12,var11); 
        Pop.Rxx{2} = Pop.Rxx{2} + int(gri * N{ij(4),ij(3)},rr1,var12,var11);
    end
    
    if includeL(5)
        %Pop.Rxy = Pop.Rxy + int(constructPxy(Q{3,5},Zxa_st,Zyo_s,gts,nx,ny,nZxa,nZyo),var12,var11,I(1,2));
        Pop.Rxy = Pop.Rxy + int(grs * subs(N{ij(3),ij(5)},var22,var21),rr1,var11,I(1,2));
    end
    if includeL(6) && ~options.sep(2)==1
        %Pop.Rxy = Pop.Rxy + int(int(constructPxy(Q{3,6},Zxa_st,Zya_st,gtt,nx,ny,nZxa,nZya),var22,var21,I(2,2)),var12,var11,I(1,2));
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(3),ij(6)},var22,var21),rr1,var11,I(1,2)),rr2,var21,I(2,2));
    elseif includeL(6) && options.sep(2)==1
        %Pop.Rxy = Pop.Rxy + int(int(constructPxy(Q{3,6},Zxa_st,Zya_st,gtt,nx,ny,nZxa,nZya),var22,I(2,1),I(2,2)),var12,var11,I(1,2));
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(3),ij(6)},var22,var21),rr1,var11,I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        %Pop.Rxy = Pop.Rxy + + int(int(constructPxy(Q{3,7},Zxa_st,Zyb_st,gtt,nx,ny,nZxa,nZyb),var22,I(2,1),var21),var12,var11,I(1,2));
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(3),ij(7)},var22,var21),rr1,var11,I(1,2)),rr2,I(2,1),var21);
    end
    
    if includeL(8)
        %Pop.Rx2{3} = Pop.Rx2{3} + constructPx2b(Q{3,8},Zxa_st,Z2oo_ss,gts,nx,n2,nZxa,nZ2oo);
        Pop.Rx2{3} = Pop.Rx2{3} + gts * subs(N{ij(3),ij(8)},[rr1;var22],[var12;var21]);
    end
    if includeL(11) && ~options.sep(4)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2b(Q{3,11},Zxa_st,Z2oa_sst,gtt,nx,n2,nZxa,nZ2oa),var22,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(gtr * subs(N{ij(3),ij(11)},[rr1;var22],[var12;var21]),rr2,var21,I(2,2));
    elseif includeL(11) && options.sep(4)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2b(Q{3,11},Zxa_st,Z2oa_sst,gtt,nx,n2,nZxa,nZ2oa),var22,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(gtr * subs(N{ij(3),ij(11)},[rr1;var22],[var12;var21]),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2b(Q{3,12},Zxa_st,Z2ob_sst,gtt,nx,n2,nZxa,nZ2ob),var22,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(gtr * subs(N{ij(3),ij(12)},[rr1;var22],[var12;var21]),rr2,I(2,1),var21);
    end
    
    if includeL(9) && ~options.sep(3)==1
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2c(Q{3,9},Zxa_st,Z2ao_sst,grs,nx,n2,nZxa,nZ2ao,rr1),rr1,var11,I(1,2));
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2c(Q{3,9},Zxa_st,Z2ao_sst,grs,nx,n2,nZxa,nZ2ao,rr1),rr1,var12,I(1,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(grs * subs(N{ij(3),ij(9)},var22,var21),rr1,var11,I(1,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(grs * subs(N{ij(3),ij(9)},var22,var21),rr1,var12,I(1,2));
    elseif includeL(9) && options.sep(3)==1
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(grs * subs(N{ij(3),ij(9)},var22,var21),rr1,var11,I(1,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(grs * subs(N{ij(3),ij(9)},var22,var21),rr1,var11,I(1,2));
    end
    if includeL(13) && ~options.sep(5)==1 && ~options.sep(6)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{3,13},Zxa_st,Z2aa_sstt,grt,nx,n2,nZxa,nZ2aa,rr1),rr1,var11,I(1,2)),var22,var21,I(2,2));
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,13},Zxa_st,Z2aa_sstt,grt,nx,n2,nZxa,nZ2aa,rr1),rr1,var12,I(1,2)),var22,var21,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var11,I(1,2)),rr2,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var12,I(1,2)),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var11,I(1,2)),rr2,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var11,I(1,2)),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{3,13},Zxa_st,Z2aa_sstt,grt,nx,n2,nZxa,nZ2aa,rr1),rr1,var11,I(1,2)),var22,I(2,1),I(2,2));
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,13},Zxa_st,Z2aa_sstt,grt,nx,n2,nZxa,nZ2aa,rr1),rr1,var12,I(1,2)),var22,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var11,I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var12,I(1,2)),rr2,I(2,1),I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var11,I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var11,I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~options.sep(5)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{3,15},Zxa_st,Z2ab_sstt,grt,nx,n2,nZxa,nZ2ab,rr1),rr1,var11,I(1,2)),var22,I(2,1),var21);
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,15},Zxa_st,Z2ab_sstt,grt,nx,n2,nZxa,nZ2ab,rr1),rr1,var12,I(1,2)),var22,I(2,1),var21);
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(15)},var22,var21),rr1,var11,I(1,2)),rr2,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(15)},var22,var21),rr1,var12,I(1,2)),rr2,I(2,1),var21);
    elseif includeL(15) && options.sep(5)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(15)},var22,var21),rr1,var11,I(1,2)),rr2,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(15)},var22,var21),rr1,var11,I(1,2)),rr2,I(2,1),var21);
    end  
    
    if includeL(10)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2c(Q{3,10},Zxa_st,Z2bo_sst,grs,nx,n2,nZxa,nZ2bo,rr1),rr1,var11,var12);
        Pop.Rx2{3} = Pop.Rx2{3} + int(grs * subs(N{ij(3),ij(10)},var22,var21),rr1,var11,var12);
    end
    if includeL(14) && ~options.sep(6)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,14},Zxa_st,Z2ba_sstt,grt,nx,n2,nZxa,nZ2ba,rr1),rr1,var11,var12),var22,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(14)},var22,var21),rr1,var11,var12),rr2,var21,I(2,2));
    elseif includeL(14) && options.sep(6)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,14},Zxa_st,Z2ba_sstt,grt,nx,n2,nZxa,nZ2ba,rr1),rr1,var11,var12),var22,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(14)},var22,var21),rr1,var11,var12),rr2,I(2,1),I(2,2));
    end
    if includeL(16)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,16},Zxa_st,Z2bb_sstt,grt,nx,n2,nZxa,nZ2bb,rr1),rr1,var11,var12),var22,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(16)},var22,var21),rr1,var11,var12),rr2,I(2,1),var21);
    end       
    
elseif includeL(3) && options.sep(1)==1
    %Pop.Rxx{2} = Pop.Rxx{2} + int(constructCA1(Q{3,3},Zxa_st,gri,nx,nZxa,var11,var12,rr1),rr1,I(1,1),I(1,2)); 
    Pop.Rxx{2} = Pop.Rxx{2} + int(gri * N{ij(3),ij(3)},rr1,I(1,1),I(1,2));
    if includeL(5)
        %Pop.Rxy = Pop.Rxy + int(constructPxy(Q{3,5},Zxa_st,Zyo_s,gts,nx,ny,nZxa,nZyo),var12,I(1,1),I(1,2));
        Pop.Rxy = Pop.Rxy + int(grs * subs(N{ij(3),ij(5)},var22,var21),rr1,I(1,1),I(1,2));
    end
    if includeL(6) && ~options.sep(2)==1
        %Pop.Rxy = Pop.Rxy + int(int(constructPxy(Q{3,6},Zxa_st,Zya_st,gtt,nx,ny,nZxa,nZya),var22,var21,I(2,2)),var12,I(1,1),I(1,2));
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(3),ij(6)},var22,var21),rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
    elseif includeL(6) && options.sep(2)==1
        %Pop.Rxy = Pop.Rxy + int(int(constructPxy(Q{3,6},Zxa_st,Zya_st,gtt,nx,ny,nZxa,nZya),var22,I(2,1),I(2,2)),var12,I(1,1),I(1,2));
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(3),ij(6)},var22,var21),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        %Pop.Rxy = Pop.Rxy + + int(int(constructPxy(Q{3,7},Zxa_st,Zyb_st,gtt,nx,ny,nZxa,nZyb),var22,I(2,1),var21),var12,I(1,1),I(1,2));
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(3),ij(7)},var22,var21),rr1,I(1,1),I(1,2)),rr2,I(2,1),var21);
    end
    
    if includeL(8)
        %Pop.Rx2{3} = Pop.Rx2{3} + constructPx2b(Q{3,8},Zxa_st,Z2oo_ss,gts,nx,n2,nZxa,nZ2oo);
        %Pop.Rx2{2} = Pop.Rx2{2} + constructPx2b(Q{3,8},Zxa_st,Z2oo_ss,gts,nx,n2,nZxb,nZ2oo);
        Pop.Rx2{3} = Pop.Rx2{3} + gts * subs(N{ij(3),ij(8)},[rr1;var22],[var12;var21]);
        Pop.Rx2{2} = Pop.Rx2{2} + gts * subs(N{ij(3),ij(8)},[rr1;var22],[var12;var21]);
    end
    if includeL(11) && ~options.sep(4)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2b(Q{3,11},Zxa_st,Z2oa_sst,gtt,nx,n2,nZxa,nZ2oa),var22,var21,I(2,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2b(Q{3,11},Zxa_st,Z2oa_sst,gtt,nx,n2,nZxa,nZ2oa),var22,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(gtr * subs(N{ij(3),ij(11)},[rr1;var22],[var12;var21]),rr2,var21,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(gtr * subs(N{ij(3),ij(11)},[rr1;var22],[var12;var21]),rr2,var21,I(2,2));
    elseif includeL(11) && options.sep(4)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2b(Q{3,11},Zxa_st,Z2oa_sst,gtt,nx,n2,nZxa,nZ2oa),var22,I(2,1),I(2,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2b(Q{3,11},Zxa_st,Z2oa_sst,gtt,nx,n2,nZxa,nZ2oa),var22,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(gtr * subs(N{ij(3),ij(11)},[rr1;var22],[var12;var21]),rr2,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(gtr * subs(N{ij(3),ij(11)},[rr1;var22],[var12;var21]),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2b(Q{3,12},Zxa_st,Z2ob_sst,gtt,nx,n2,nZxa,nZ2ob),var22,I(2,1),var21);
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2b(Q{3,12},Zxa_st,Z2ob_sst,gtt,nx,n2,nZxa,nZ2ob),var22,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(gtr * subs(N{ij(3),ij(12)},[rr1;var22],[var12;var21]),rr2,I(2,1),var21);
        Pop.Rx2{2} = Pop.Rx2{2} + int(gtr * subs(N{ij(3),ij(12)},[rr1;var22],[var12;var21]),rr2,I(2,1),var21);
    end
    
    if includeL(9) && ~options.sep(3)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2c(Q{3,9},Zxa_st,Z2ao_sst,grs,nx,n2,nZxa,nZ2ao,rr1),rr1,var11,I(1,2));
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2c(Q{3,9},Zxa_st,Z2ao_sst,grs,nx,n2,nZxa,nZ2ao,rr1),rr1,var12,I(1,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2c(Q{4,9},Zxb_st,Z2ao_sst,grs,nx,n2,nZxb,nZ2ao,rr1),rr1,var12,var11);
        Pop.Rx2{2} = Pop.Rx2{2} + int(grs * subs(N{ij(3),ij(9)},var22,var21),rr1,var12,I(1,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(grs * subs(N{ij(3),ij(9)},var22,var21),rr1,var12,I(1,2));
    elseif options.sep(3)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(grs * subs(N{ij(3),ij(9)},var22,var21),rr1,I(1,1),I(1,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(grs * subs(N{ij(3),ij(9)},var22,var21),rr1,I(1,1),I(1,2));
    end
    if includeL(13) && ~options.sep(5)==1 && ~options.sep(6)==1
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{3,13},Zxa_st,Z2aa_sstt,grt,nx,n2,nZxa,nZ2aa,rr1),rr1,var11,I(1,2)),var22,var21,I(2,2));
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,13},Zxa_st,Z2aa_sstt,grt,nx,n2,nZxa,nZ2aa,rr1),rr1,var12,I(1,2)),var22,var21,I(2,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,13},Zxb_st,Z2aa_sstt,grt,nx,n2,nZxb,nZ2aa,rr1),rr1,var12,var11),var22,var21,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var12,I(1,2)),rr2,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var12,I(1,2)),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{3,13},Zxa_st,Z2aa_sstt,grt,nx,n2,nZxa,nZ2aa,rr1),rr1,var11,I(1,2)),var22,I(2,1),I(2,2));
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,13},Zxa_st,Z2aa_sstt,grt,nx,n2,nZxa,nZ2aa,rr1),rr1,var12,I(1,2)),var22,I(2,1),I(2,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,13},Zxb_st,Z2aa_sstt,grt,nx,n2,nZxb,nZ2aa,rr1),rr1,var12,var11),var22,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var12,I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,var12,I(1,2)),rr2,I(2,1),I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(13)},var22,var21),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));        
    end
    if includeL(15) && ~options.sep(5)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{3,15},Zxa_st,Z2ab_sstt,grt,nx,n2,nZxa,nZ2ab,rr1),rr1,var11,I(1,2)),var22,I(2,1),var21);
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,15},Zxa_st,Z2ab_sstt,grt,nx,n2,nZxa,nZ2ab,rr1),rr1,var12,I(1,2)),var22,I(2,1),var21);
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,15},Zxb_st,Z2ab_sstt,grt,nx,n2,nZxb,nZ2ab,rr1),rr1,var12,var11),var22,I(2,1),var21);
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(15)},var22,var21),rr1,var12,I(1,2)),rr2,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(15)},var22,var21),rr1,var12,I(1,2)),rr2,I(2,1),var21);
    elseif includeL(15) && options.sep(5)
        %
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(15)},var22,var21),rr1,I(1,1),I(1,2)),rr2,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(15)},var22,var21),rr1,I(1,1),I(1,2)),rr2,I(2,1),var21);
    end  
    
    if includeL(10)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2c(Q{3,10},Zxa_st,Z2bo_sst,grs,nx,n2,nZxa,nZ2bo,rr1),rr1,var11,var12);
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2c(Q{4,10},Zxb_st,Z2bo_sst,grs,nx,n2,nZxb,nZ2bo,rr1),rr1,I(1,1),var11);
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2c(Q{4,10},Zxb_st,Z2bo_sst,grs,nx,n2,nZxb,nZ2bo,rr1),rr1,I(1,1),var12);
        Pop.Rx2{2} = Pop.Rx2{2} + int(grs * subs(N{ij(3),ij(10)},var22,var21),rr1,I(1,1),var12);
        Pop.Rx2{3} = Pop.Rx2{3} + int(grs * subs(N{ij(3),ij(10)},var22,var21),rr1,I(1,1),var12);
    end
    if includeL(14) && ~options.sep(6)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,14},Zxa_st,Z2ba_sstt,grt,nx,n2,nZxa,nZ2ba,rr1),rr1,var11,var12),var22,var21,I(2,2));
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{4,14},Zxb_st,Z2ba_sstt,grt,nx,n2,nZxb,nZ2ba,rr1),rr1,I(1,1),var11),var22,var21,I(2,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,14},Zxb_st,Z2ba_sstt,grt,nx,n2,nZxb,nZ2ba,rr1),rr1,I(1,1),var12),var22,var21,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(14)},var22,var21),rr1,I(1,1),var12),rr2,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(14)},var22,var21),rr1,I(1,1),var12),rr2,var21,I(2,2));
    elseif includeL(14) && options.sep(6)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,14},Zxa_st,Z2ba_sstt,grt,nx,n2,nZxa,nZ2ba,rr1),rr1,var11,var12),var22,I(2,1),I(2,2));
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{4,14},Zxb_st,Z2ba_sstt,grt,nx,n2,nZxb,nZ2ba,rr1),rr1,I(1,1),var11),var22,I(2,1),I(2,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,14},Zxb_st,Z2ba_sstt,grt,nx,n2,nZxb,nZ2ba,rr1),rr1,I(1,1),var12),var22,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(14)},var22,var21),rr1,I(1,1),var12),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(14)},var22,var21),rr1,I(1,1),var12),rr2,I(2,1),I(2,2));
    end
    if includeL(16)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{3,16},Zxa_st,Z2bb_sstt,grt,nx,n2,nZxa,nZ2bb,rr1),rr1,var11,var12),var22,I(2,1),var21);
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{4,16},Zxb_st,Z2bb_sstt,grt,nx,n2,nZxb,nZ2bb,rr1),rr1,I(1,1),var11),var22,I(2,1),var21);
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,16},Zxb_st,Z2bb_sstt,grt,nx,n2,nZxb,nZ2bb,rr1),rr1,I(1,1),var12),var22,I(2,1),var21);
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(3),ij(16)},var22,var21),rr1,I(1,1),var12),rr2,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(3),ij(16)},var22,var21),rr1,I(1,1),var12),rr2,I(2,1),var21);
    end  
    
end
if includeL(4)
    %Pop.Rxx{2} = Pop.Rxx{2} + int(constructCA1(Q{4,4},Zxa_st,gri,nx,nZxa,var11,var12,rr1),rr1,I(1,1),var12);
    Pop.Rxx{2} = Pop.Rxx{2} + int(gri * N{ij(4),ij(4)},rr1,I(1,1),var12);
    if includeL(5)
        %Pop.Rxy = Pop.Rxy + int(constructPxy(Q{4,5},Zxb_st,Zyo_s,gts,nx,ny,nZxb,nZyo),var12,I(1,1),var11);
        Pop.Rxy = Pop.Rxy + int(grs * subs(N{ij(4),ij(5)},var22,var21),rr1,I(1,1),var11);
    end
    if includeL(6) && ~options.sep(2)==1
        %Pop.Rxy = Pop.Rxy + int(int(constructPxy(Q{4,6},Zxb_st,Zya_st,gtt,nx,ny,nZxb,nZya),var22,var21,I(2,2)),var12,I(1,1),var11);
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(4),ij(6)},var22,var21),rr1,I(1,1),var11),rr2,var21,I(2,2));
    elseif includeL(6) && options.sep(2)==1
        %Pop.Rxy = Pop.Rxy + int(int(constructPxy(Q{4,6},Zxb_st,Zya_st,gtt,nx,ny,nZxb,nZya),var22,I(2,1),I(2,2)),var12,I(1,1),var11);
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(4),ij(6)},var22,var21),rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
    end
    if includeL(7)
        %Pop.Rxy = Pop.Rxy + int(int(constructPxy(Q{4,7},Zxb_st,Zyb_st,gtt,nx,ny,nZxb,nZyb),var22,I(2,1),var21),var12,I(1,1),var11);
        Pop.Rxy = Pop.Rxy + int(int(grr * subs(N{ij(4),ij(7)},var22,var21),rr1,I(1,1),var11),rr2,I(2,1),var21);
    end
    
    if includeL(8)
        %Pop.Rx2{2} = Pop.Rx2{2} + constructPx2b(Q{4,8},Zxb_st,Z2oo_ss,gts,nx,n2,nZxb,nZ2oo);
        Pop.Rx2{2} = Pop.Rx2{2} + gts * subs(N{ij(4),ij(8)},[rr1;var22],[var12;var21]);
    end
    if includeL(11) && ~options.sep(4)==1
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2b(Q{4,11},Zxb_st,Z2oa_sst,gtt,nx,n2,nZxb,nZ2oa),var22,var21,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(gtr * subs(N{ij(4),ij(11)},[rr1;var22],[var12;var21]),rr2,var21,I(2,2));
    elseif includeL(11) && options.sep(4)==1
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2b(Q{4,11},Zxb_st,Z2oa_sst,gtt,nx,n2,nZxb,nZ2oa),var22,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(gtr * subs(N{ij(4),ij(11)},[rr1;var22],[var12;var21]),rr2,I(2,1),I(2,2));
    end
    if includeL(12)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2b(Q{4,12},Zxb_st,Z2ob_sst,gtt,nx,n2,nZxb,nZ2ob),var22,I(2,1),var21);
        Pop.Rx2{2} = Pop.Rx2{2} + int(gtr * subs(N{ij(4),ij(12)},[rr1;var22],[var12;var21]),rr2,I(2,1),var21);
    end
    
    if includeL(9) && ~options.sep(3)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2c(Q{4,9},Zxb_st,Z2ao_sst,grs,nx,n2,nZxb,nZ2ao,rr1),rr1,var12,var11);
        Pop.Rx2{2} = Pop.Rx2{2} + int(grs * subs(N{ij(4),ij(9)},var22,var21),rr1,var12,var11);
    elseif includeL(9) && options.sep(3)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(grs * subs(N{ij(4),ij(9)},var22,var21),rr1,I(1,1),var11);
        Pop.Rx2{3} = Pop.Rx2{3} + int(grs * subs(N{ij(4),ij(9)},var22,var21),rr1,I(1,1),var11);
    end
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,13},Zxb_st,Z2aa_sstt,grt,nx,n2,nZxb,nZ2aa,rr1),rr1,var12,var11),var22,var21,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(13)},var22,var21),rr1,var12,var11),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(13)},var22,var21),rr1,I(1,1),var11),rr2,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(4),ij(13)},var22,var21),rr1,I(1,1),var11),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(13)},var22,var21),rr1,var12,var11),rr2,I(2,1),I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,13},Zxb_st,Z2aa_sstt,grt,nx,n2,nZxb,nZ2aa,rr1),rr1,var12,var11),var22,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(13)},var22,var21),rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(4),ij(13)},var22,var21),rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~options.sep(5)
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,15},Zxb_st,Z2ab_sstt,grt,nx,n2,nZxb,nZ2ab,rr1),rr1,var12,var11),var22,I(2,1),var21);
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(15)},var22,var21),rr1,var12,var11),rr2,I(2,1),var21);
    elseif includeL(15) && options.sep(5)
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(15)},var22,var21),rr1,I(1,1),var11),rr2,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(4),ij(15)},var22,var21),rr1,I(1,1),var11),rr2,I(2,1),var21);
    end     
    
    if includeL(10)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(constructPx2c(Q{4,10},Zxb_st,Z2bo_sst,grs,nx,n2,nZxb,nZ2bo,rr1),rr1,I(1,1),var11);
        %Pop.Rx2{2} = Pop.Rx2{2} + int(constructPx2c(Q{4,10},Zxb_st,Z2bo_sst,grs,nx,n2,nZxb,nZ2bo,rr1),rr1,I(1,1),var12);
        Pop.Rx2{2} = Pop.Rx2{2} + int(grs * subs(N{ij(4),ij(10)},var22,var21),rr1,I(1,1),var12);
        Pop.Rx2{3} = Pop.Rx2{3} + int(grs * subs(N{ij(4),ij(10)},var22,var21),rr1,I(1,1),var11);
    end
    if includeL(14) && ~options.sep(6)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{4,14},Zxb_st,Z2ba_sstt,grt,nx,n2,nZxb,nZ2ba,rr1),rr1,I(1,1),var11),var22,var21,I(2,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,14},Zxb_st,Z2ba_sstt,grt,nx,n2,nZxb,nZ2ba,rr1),rr1,I(1,1),var12),var22,var21,I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(14)},var22,var21),rr1,I(1,1),var12),rr2,var21,I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(4),ij(14)},var22,var21),rr1,I(1,1),var11),rr2,var21,I(2,2));
    elseif includeL(14) && options.sep(6)==1
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{4,14},Zxb_st,Z2ba_sstt,grt,nx,n2,nZxb,nZ2ba,rr1),rr1,I(1,1),var11),var22,I(2,1),I(2,2));
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,14},Zxb_st,Z2ba_sstt,grt,nx,n2,nZxb,nZ2ba,rr1),rr1,I(1,1),var12),var22,I(2,1),I(2,2));
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(14)},var22,var21),rr1,I(1,1),var12),rr2,I(2,1),I(2,2));
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(4),ij(14)},var22,var21),rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
    end
    if includeL(16)
        %Pop.Rx2{3} = Pop.Rx2{3} + int(int(constructPx2c(Q{4,16},Zxb_st,Z2bb_sstt,grt,nx,n2,nZxb,nZ2bb,rr1),rr1,I(1,1),var11),var22,I(2,1),var21);
        %Pop.Rx2{2} = Pop.Rx2{2} + int(int(constructPx2c(Q{4,16},Zxb_st,Z2bb_sstt,grt,nx,n2,nZxb,nZ2bb,rr1),rr1,I(1,1),var12),var22,I(2,1),var21);
        Pop.Rx2{2} = Pop.Rx2{2} + int(int(grr * subs(N{ij(4),ij(16)},var22,var21),rr1,I(1,1),var12),rr2,I(2,1),var21);
        Pop.Rx2{3} = Pop.Rx2{3} + int(int(grr * subs(N{ij(4),ij(16)},var22,var21),rr1,I(1,1),var11),rr2,I(2,1),var21);
    end        
    
end

if options.sep(1)==1
    Pop.Rxx{3} = Pop.Rxx{2};
else
    Pop.Rxx{3} = var_swap(Pop.Rxx{2},var11,var12).';
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Build Ryy and Ry2
if includeL(5)
    %Pop.Ryy{1} = constructR0(Q{5,5},Zyo_s,gis,ny,nZyo); 
    Pop.Ryy{1} = gis * subs(N{ij(5),ij(5)},var22,var21); 
    if includeL(6) && ~options.sep(2)==1
        %Pop.Ryy{2} = Pop.Ryy{2} + constructA1(Q{5,6},Zyo_s,Zya_st,gis,ny,nZyo,nZya,var21,var22); 
        Pop.Ryy{2} = Pop.Ryy{2} + gis * subs(N{ij(5),ij(6)},rr2,var21); 
    elseif includeL(6) && options.sep(1)==1
        %Pop.Ryy{2} = Pop.Ryy{2} + constructA1(Q{5,6},Zyo_s,Zya_st,gis,ny,nZyo,nZya,var21,var22) + constructB1(Q{6,5},Zyo_s,Zya_st,git,ny,nZyo,nZya,var21,var22);
        Pop.Ryy{2} = Pop.Ryy{2} + gis * subs(N{ij(5),ij(6)},rr2,var21) + git * subs(N{ij(6),ij(5)},rr2,var22); 
    end
    if includeL(7)
        %Pop.Ryy{2} = Pop.Ryy{2} + constructB1(Q{7,5},Zyo_s,Zya_st,git,ny,nZyo,nZya,var21,var22);
        Pop.Ryy{2} = Pop.Ryy{2} + git * subs(N{ij(7),ij(5)},rr2,var22); 
    end
    
    if includeL(8)
        %Pop.Ry2{1} = Pop.Ry2{1} + constructPy2a(Q{5,8},Zyo_s,Z2oo_ss,gss,ny,n2,nZyo,nZ2oo,[0;0]);
        Pop.Ry2{1} = Pop.Ry2{1} + gss * subs(N{ij(5),ij(8)},var2,var1); 
    end
    if includeL(9) && ~options.sep(3)==1
        %Pop.Ry2{1} = Pop.Ry2{1} + int(constructPy2a(Q{5,9},Zyo_s,Z2ao_sst,gts,ny,n2,nZyo,nZ2ao,[1;0]),var12,var11,I(1,2));
        Pop.Ry2{1} = Pop.Ry2{1} + int(grs * subs(N{ij(5),ij(9)},var2,var1),rr1,var11,I(1,2)); 
    elseif includeL(9) && options.sep(3)==1
        %Pop.Ry2{1} = Pop.Ry2{1} + int(constructPy2a(Q{5,9},Zyo_s,Z2ao_sst,gts,ny,n2,nZyo,nZ2ao,[1;0]),var12,I(1,1),I(1,2));
        Pop.Ry2{1} = Pop.Ry2{1} + int(grs * subs(N{ij(5),ij(9)},var2,var1),rr1,I(1,1),I(1,2)); 
    end
    if includeL(10)
        %Pop.Ry2{1} = Pop.Ry2{1} + int(constructPy2a(Q{5,10},Zyo_s,Z2bo_sst,gts,ny,n2,nZyo,nZ2bo,[1;0]),var12,I(1,1),var11);
        Pop.Ry2{1} = Pop.Ry2{1} + int(grs * subs(N{ij(5),ij(10)},var2,var1),rr1,I(1,1),var11); 
    end
    
    if includeL(11) && ~options.sep(4)
        %Pop.Ry2{2} = Pop.Ry2{2} + constructPy2a(Q{5,11},Zyo_s,Z2oa_sst,gss,ny,n2,nZyo,nZ2oa,[0;0]);
        Pop.Ry2{2} = Pop.Ry2{2} + gss * subs(N{ij(5),ij(11)},[var12;rr2],var1);
    elseif includeL(11) && options.sep(4)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + gss * subs(N{ij(5),ij(11)},[var12;rr2],var1);
        Pop.Ry2{3} = Pop.Ry2{3} + gss * subs(N{ij(5),ij(11)},[var12;rr2],var1);
    end
    if includeL(13) && ~options.sep(5)==1 && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2a(Q{5,13},Zyo_s,Z2aa_sstt,gts,ny,n2,nZyo,nZ2aa,[1;0]),var12,var11,I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(grs * subs(N{ij(5),ij(13)},[var12;rr2],var1),rr1,var11,I(1,2));
    elseif includeL(13) && ~options.sep(6)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(grs * subs(N{ij(5),ij(13)},[var12;rr2],var1),rr1,I(1,1),I(1,2));
    elseif includeL(13) && ~options.sep(5)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(grs * subs(N{ij(5),ij(13)},[var12;rr2],var1),rr1,var11,I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(grs * subs(N{ij(5),ij(13)},[var12;rr2],var1),rr1,var11,I(1,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2a(Q{5,13},Zyo_s,Z2aa_sstt,gts,ny,n2,nZyo,nZ2aa,[1;0]),var12,I(1,1),I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(grs * subs(N{ij(5),ij(13)},[var12;rr2],var1),rr1,I(1,1),I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(grs * subs(N{ij(5),ij(13)},[var12;rr2],var1),rr1,I(1,1),I(1,2));
    end
    if includeL(14) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2a(Q{5,14},Zyo_s,Z2ba_sstt,gts,ny,n2,nZyo,nZ2ba,[1;0]),var12,I(1,1),var11);
        Pop.Ry2{2} = Pop.Ry2{2} + int(grs * subs(N{ij(5),ij(14)},[var12;rr2],var1),rr1,I(1,1),var11);
    elseif includeL(14) && options.sep(6)
        Pop.Ry2{2} = Pop.Ry2{2} + int(grs * subs(N{ij(5),ij(14)},[var12;rr2],var1),rr1,I(1,1),var11);
        Pop.Ry2{3} = Pop.Ry2{3} + int(grs * subs(N{ij(5),ij(14)},[var12;rr2],var1),rr1,I(1,1),var11);
    end
    
    if includeL(12)
        %Pop.Ry2{3} = Pop.Ry2{3} + constructPy2a(Q{5,12},Zyo_s,Z2ob_sst,gss,ny,n2,nZyo,nZ2ob,[0;0]);
        Pop.Ry2{3} = Pop.Ry2{3} + gss * subs(N{ij(5),ij(12)},[var12;rr2],var1);
    end
    if includeL(15) && ~options.sep(5)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2a(Q{5,15},Zyo_s,Z2ab_sstt,gts,ny,n2,nZyo,nZ2ab,[1;0]),var12,var11,I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(grs * subs(N{ij(5),ij(15)},[var12;rr2],var1),rr1,var11,I(1,2));
    elseif includeL(15) && options.sep(5)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2a(Q{5,15},Zyo_s,Z2ab_sstt,gts,ny,n2,nZyo,nZ2ab,[1;0]),var12,I(1,1),I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(grs * subs(N{ij(5),ij(15)},[var12;rr2],var1),rr1,I(1,1),I(1,2));
    end
    if includeL(16)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2a(Q{5,16},Zyo_s,Z2bb_sstt,gts,ny,n2,nZyo,nZ2bb,[1;0]),var12,I(1,1),var11);
        Pop.Ry2{3} = Pop.Ry2{3} + int(grs * subs(N{ij(5),ij(16)},[var12;rr2],var1),rr1,I(1,1),var11);
    end 
end

if includeL(6) && ~options.sep(2)==1
    %Pop.Ryy{2} = Pop.Ryy{2} + int(constructCA1(Q{6,6},Zya_st,gir,ny,nZya,var21,var22,rr2),rr2,var21,I(2,2));
    Pop.Ryy{2} = Pop.Ryy{2} + int(gir * N{ij(6),ij(6)},rr2,var21,I(2,2));
    if includeL(7)
        %Pop.Ryy{2} = Pop.Ryy{2} + int(constructCA1(Q{7,6},Zya_st,gir,ny,nZya,var21,var22,rr2),rr2,var22,var21);
        Pop.Ryy{2} = Pop.Ryy{2} + int(gir * N{ij(7),ij(6)},rr2,var22,var21);
    end
    
    if includeL(8)
        %Pop.Ry2{3} = Pop.Ry2{3} + constructPy2b(Q{6,8},Zya_st,Z2oo_ss,gst,ny,n2,nZya,nZ2oo);
        Pop.Ry2{3} = Pop.Ry2{3} + gst * subs(N{ij(6),ij(8)},[var12;rr2],[var11;var22]);
    end
    if includeL(9) && ~options.sep(3)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2b(Q{6,9},Zya_st,Z2ao_sst,gtt,ny,n2,nZya,nZ2ao),var12,var11,I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(grt * subs(N{ij(6),ij(9)},[var12;rr2],[var11;var22]),rr1,var11,I(1,2));
    elseif includeL(9) && options.sep(3)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2b(Q{6,9},Zya_st,Z2ao_sst,gtt,ny,n2,nZya,nZ2ao),var12,I(1,1),I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(grt * subs(N{ij(6),ij(9)},[var12;rr2],[var11;var22]),rr1,I(1,1),I(1,2));
    end
    if includeL(10)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2b(Q{6,10},Zya_st,Z2bo_sst,gtt,ny,n2,nZya,nZ2bo),var12,I(1,1),var11);
        Pop.Ry2{3} = Pop.Ry2{3} + int(grt * subs(N{ij(6),ij(10)},[var12;rr2],[var11;var22]),rr1,I(1,1),var11);
    end
    
    if includeL(11) && ~options.sep(4)==1
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2c(Q{6,11},Zya_st,Z2oa_sst,gsr,ny,n2,nZya,nZ2oa,rr2),rr2,var21,I(2,2));
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2c(Q{6,11},Zya_st,Z2oa_sst,gsr,ny,n2,nZya,nZ2oa,rr2),rr2,var22,I(2,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(gsr * subs(N{ij(6),ij(11)},var12,var11),rr2,var21,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(gsr * subs(N{ij(6),ij(11)},var12,var11),rr2,var22,I(2,2));
    elseif includeL(11) && options.sep(4)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2c(Q{6,11},Zya_st,Z2oa_sst,gsr,ny,n2,nZya,nZ2oa,rr2),rr2,var21,I(2,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(gsr * subs(N{ij(6),ij(11)},var12,var11),rr2,var21,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(gsr * subs(N{ij(6),ij(11)},var12,var11),rr2,var21,I(2,2));
    end   
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{6,13},Zya_st,Z2aa_sstt,gtr,ny,n2,nZya,nZ2aa,rr2),rr2,var21,I(2,2)),var12,var11,I(1,2));
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,13},Zya_st,Z2aa_sstt,gtr,ny,n2,nZya,nZ2aa,rr2),rr2,var22,I(2,2)),var12,var11,I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,var21,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,var22,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{6,13},Zya_st,Z2aa_sstt,gtr,ny,n2,nZya,nZ2aa,rr2),rr2,var21,I(2,2)),var12,I(1,1),I(1,2));
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,13},Zya_st,Z2aa_sstt,gtr,ny,n2,nZya,nZ2aa,rr2),rr2,var22,I(2,2)),var12,I(1,1),I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,var22,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,var21,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,var21,I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
    end
    if includeL(14) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{6,14},Zya_st,Z2ba_sstt,gtr,ny,n2,nZya,nZ2ba,rr2),rr2,var21,I(2,2)),var12,I(1,1),var11);
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,14},Zya_st,Z2ba_sstt,gtr,ny,n2,nZya,nZ2ba,rr2),rr2,var22,I(2,2)),var12,I(1,1),var11);
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,var21,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,var22,I(2,2));
    elseif includeL(14) && options.sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,var21,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,var21,I(2,2));
    end
    
    if includeL(12)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2c(Q{6,12},Zya_st,Z2ob_sst,gsr,ny,n2,nZya,nZ2ob,rr2),rr2,var21,var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(gsr * subs(N{ij(6),ij(12)},var12,var11),rr2,var21,var22);
    end
    if includeL(15) && ~options.sep(5)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,15},Zya_st,Z2ab_sstt,gtr,ny,n2,nZya,nZ2ab,rr2),rr2,var21,var22),var12,var11,I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(15)},var12,var11),rr1,var11,I(1,2)),rr2,var21,var22);
    elseif includeL(15) && options.sep(5)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,15},Zya_st,Z2ab_sstt,gtr,ny,n2,nZya,nZ2ab,rr2),rr2,var21,var22),var12,I(1,1),I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(15)},var12,var11),rr1,I(1,1),I(1,2)),rr2,var21,var22);
    end
    if includeL(16)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,16},Zya_st,Z2bb_sstt,gtr,ny,n2,nZya,nZ2bb,rr2),rr2,var21,var22),var12,I(1,1),var11);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(16)},var12,var11),rr1,I(1,1),var11),rr2,var21,var22);
    end       
   
elseif includeL(6) && options.sep(2)==1
    %Pop.Ryy{2} = Pop.Ryy{2} + int(constructCA1(Q{6,6},Zya_st,gir,ny,nZya,var21,var22,rr2),rr2,I(2,1),I(2,2)); 
    Pop.Ryy{2} = Pop.Ryy{2} + int(gir * N{ij(6),ij(6)},rr2,I(2,1),I(2,2));
    
    if includeL(8)
        %Pop.Ry2{3} = Pop.Ry2{3} + constructPy2b(Q{6,8},Zya_st,Z2oo_ss,gst,ny,n2,nZya,nZ2oo);
        %Pop.Ry2{2} = Pop.Ry2{2} + constructPy2b(Q{6,8},Zya_st,Z2oo_ss,gst,ny,n2,nZya,nZ2oo);
        Pop.Ry2{2} = Pop.Ry2{2} + gst * subs(N{ij(6),ij(8)},[var12;rr2],[var11;var22]);
        Pop.Ry2{3} = Pop.Ry2{3} + gst * subs(N{ij(6),ij(8)},[var12;rr2],[var11;var22]);
    end
    if includeL(9) && ~options.sep(3)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2b(Q{6,9},Zya_st,Z2ao_sst,gtt,ny,n2,nZya,nZ2ao),var12,var11,I(1,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2b(Q{6,9},Zya_st,Z2ao_sst,gtt,ny,n2,nZya,nZ2ao),var12,var11,I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(grt * subs(N{ij(6),ij(9)},[var12;rr2],[var11;var22]),rr1,var11,I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(grt * subs(N{ij(6),ij(9)},[var12;rr2],[var11;var22]),rr1,var11,I(1,2));
    elseif includeL(9) && options.sep(3)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2b(Q{6,9},Zya_st,Z2ao_sst,gtt,ny,n2,nZya,nZ2ao),var12,I(1,1),I(1,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2b(Q{6,9},Zya_st,Z2ao_sst,gtt,ny,n2,nZya,nZ2ao),var12,I(1,1),I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(grt * subs(N{ij(6),ij(9)},[var12;rr2],[var11;var22]),rr1,I(1,1),I(1,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(grt * subs(N{ij(6),ij(9)},[var12;rr2],[var11;var22]),rr1,I(1,1),I(1,2));
    end
    if includeL(10)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2b(Q{6,10},Zya_st,Z2bo_sst,gtt,ny,n2,nZya,nZ2bo),var12,I(1,1),var11);
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2b(Q{6,10},Zya_st,Z2bo_sst,gtt,ny,n2,nZya,nZ2bo),var12,I(1,1),var11);
        Pop.Ry2{2} = Pop.Ry2{2} + int(grt * subs(N{ij(6),ij(10)},[var12;rr2],[var11;var22]),rr1,I(1,1),var11);
        Pop.Ry2{3} = Pop.Ry2{3} + int(grt * subs(N{ij(6),ij(10)},[var12;rr2],[var11;var22]),rr1,I(1,1),var11);
    end
    
    if includeL(11) && ~options.sep(4)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2c(Q{6,11},Zya_st,Z2oa_sst,gsr,ny,n2,nZya,nZ2oa,rr2),rr2,var21,I(2,2));
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2c(Q{6,11},Zya_st,Z2oa_sst,gsr,ny,n2,nZya,nZ2oa,rr2),rr2,var22,I(2,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2c(Q{7,11},Zyb_st,Z2oa_sst,gsr,ny,n2,nZyb,nZ2oa,rr2),rr2,var22,var21);
        Pop.Ry2{2} = Pop.Ry2{2} + int(gsr * subs(N{ij(6),ij(11)},var12,var11),rr2,var22,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(gsr * subs(N{ij(6),ij(11)},var12,var11),rr2,var22,I(2,2));
    elseif includeL(11) && options.sep(4)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(gsr * subs(N{ij(6),ij(11)},var12,var11),rr2,I(2,1),I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(gsr * subs(N{ij(6),ij(11)},var12,var11),rr2,I(2,1),I(2,2));
    end  
    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{6,13},Zya_st,Z2aa_sstt,gtr,ny,n2,nZya,nZ2aa,rr2),rr2,var21,I(2,2)),var12,var11,I(1,2));
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,13},Zya_st,Z2aa_sstt,gtr,ny,n2,nZya,nZ2aa,rr2),rr2,var22,I(2,2)),var12,var11,I(1,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,13},Zyb_st,Z2aa_sstt,gtr,ny,n2,nZyb,nZ2aa,rr2),rr2,var22,var21),var12,var11,I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,var22,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,var22,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{6,13},Zya_st,Z2aa_sstt,gtr,ny,n2,nZya,nZ2aa,rr2),rr2,var21,I(2,2)),var12,I(1,1),I(1,2));
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,13},Zya_st,Z2aa_sstt,gtr,ny,n2,nZya,nZ2aa,rr2),rr2,var22,I(2,2)),var12,I(1,1),I(1,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,13},Zyb_st,Z2aa_sstt,gtr,ny,n2,nZyb,nZ2aa,rr2),rr2,var22,var21),var12,I(1,1),I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,var22,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,var22,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,I(2,1),I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    end
    if includeL(14) && ~options.sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,var22,I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,var22,I(2,2));
    elseif includeL(14) && options.sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
    end
    
    if includeL(12)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2c(Q{6,12},Zya_st,Z2ob_sst,gsr,ny,n2,nZya,nZ2ob,rr2),rr2,var21,var22);
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2c(Q{7,12},Zyb_st,Z2ob_sst,gsr,ny,n2,nZyb,nZ2ob,rr2),rr2,I(2,1),var21);
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2c(Q{7,12},Zyb_st,Z2ob_sst,gsr,ny,n2,nZyb,nZ2ob,rr2),rr2,I(2,1),var22);
        Pop.Ry2{2} = Pop.Ry2{2} + int(gsr * subs(N{ij(6),ij(12)},var12,var11),rr2,I(2,1),var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(gsr * subs(N{ij(6),ij(12)},var12,var11),rr2,I(2,1),var22);
    end
    if includeL(15) && ~options.sep(5)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,15},Zya_st,Z2ab_sstt,gtr,ny,n2,nZya,nZ2ab,rr2),rr2,var21,var22),var12,var11,I(1,2));
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{7,15},Zyb_st,Z2ab_sstt,gtr,ny,n2,nZyb,nZ2ab,rr2),rr2,I(2,1),var21),var12,var11,I(1,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,15},Zyb_st,Z2ab_sstt,gtr,ny,n2,nZyb,nZ2ab,rr2),rr2,I(2,1),var22),var12,var11,I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(15)},var12,var11),rr1,var11,I(1,2)),rr2,I(2,1),var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(15)},var12,var11),rr1,var11,I(1,2)),rr2,I(2,1),var22);
    elseif includeL(15) && options.sep(5)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,15},Zya_st,Z2ab_sstt,gtr,ny,n2,nZya,nZ2ab,rr2),rr2,var21,var22),var12,I(1,1),I(1,2));
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{7,15},Zyb_st,Z2ab_sstt,gtr,ny,n2,nZyb,nZ2ab,rr2),rr2,I(2,1),var21),var12,I(1,1),I(1,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,15},Zyb_st,Z2ab_sstt,gtr,ny,n2,nZyb,nZ2ab,rr2),rr2,I(2,1),var22),var12,I(1,1),I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(15)},var12,var11),rr1,I(1,1),I(1,2)),rr2,I(2,1),var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(15)},var12,var11),rr1,I(1,1),I(1,2)),rr2,I(2,1),var22);
    end
    if includeL(16)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{6,16},Zya_st,Z2bb_sstt,gtr,ny,n2,nZya,nZ2bb,rr2),rr2,var21,var22),var12,I(1,1),var11);
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{7,16},Zyb_st,Z2bb_sstt,gtr,ny,n2,nZyb,nZ2bb,rr2),rr2,I(2,1),var21),var12,I(1,1),var11);
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,16},Zyb_st,Z2bb_sstt,gtr,ny,n2,nZyb,nZ2bb,rr2),rr2,I(2,1),var22),var12,I(1,1),var11);
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(6),ij(16)},var12,var11),rr1,I(1,1),var11),rr2,I(2,1),var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(6),ij(16)},var12,var11),rr1,I(1,1),var11),rr2,I(2,1),var22);
    end   
end

if includeL(7)
    %Pop.Ryy{2} = Pop.Ryy{2} + int(constructCA1(Q{7,7},Zya_st,gir,ny,nZya,var21,var22,rr2),rr2,I(2,1),var22);
    Pop.Ryy{2} = Pop.Ryy{2} + int(gir * N{ij(7),ij(7)},rr2,I(2,1),var22);
    
    if includeL(8)
        %Pop.Ry2{2} = Pop.Ry2{2} + constructPy2b(Q{7,8},Zyb_st,Z2oo_ss,gst,ny,n2,nZyb,nZ2oo);
        Pop.Ry2{2} = Pop.Ry2{2} + gst * subs(N{ij(7),ij(8)},[var12;rr2],[var11;var22]);
    end
    if includeL(9) && ~options.sep(3)==1
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2b(Q{7,9},Zyb_st,Z2ao_sst,gtt,ny,n2,nZyb,nZ2ao),var12,var11,I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(grt * subs(N{ij(7),ij(9)},[var12;rr2],[var11;var22]),rr1,var11,I(1,2));
    elseif includeL(9) && options.sep(3)==1
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2b(Q{7,9},Zyb_st,Z2ao_sst,gtt,ny,n2,nZyb,nZ2ao),var12,I(1,1),I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(grt * subs(N{ij(7),ij(9)},[var12;rr2],[var11;var22]),rr1,I(1,1),I(1,2));
    end
    if includeL(10)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2b(Q{7,10},Zyb_st,Z2bo_sst,gtt,ny,n2,nZyb,nZ2bo),var12,I(1,1),var11);
        Pop.Ry2{2} = Pop.Ry2{2} + int(grt * subs(N{ij(7),ij(10)},[var12;rr2],[var11;var22]),rr1,I(1,1),var11);
    end
    
    if includeL(11) && ~options.sep(4)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2c(Q{7,11},Zyb_st,Z2oa_sst,gsr,ny,n2,nZyb,nZ2oa,rr2),rr2,var22,var21);
        Pop.Ry2{2} = Pop.Ry2{2} + int(gsr * subs(N{ij(7),ij(11)},var12,var11),rr2,var22,var21);
    elseif includeL(11) && options.sep(4)
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(gsr * subs(N{ij(7),ij(11)},var12,var11),rr2,I(1,1),var21);
        Pop.Ry2{3} = Pop.Ry2{3} + int(gsr * subs(N{ij(7),ij(11)},var12,var11),rr2,I(1,1),var21);
    end    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,13},Zyb_st,Z2aa_sstt,gtr,ny,n2,nZyb,nZ2aa,rr2),rr2,var22,var21),var12,var11,I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,var22,var21);
    elseif includeL(13) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,13},Zyb_st,Z2aa_sstt,gtr,ny,n2,nZyb,nZ2aa,rr2),rr2,var22,var21),var12,I(1,1),I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,var22,var21);
    elseif includeL(13) && ~options.sep(5)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,I(2,1),var21);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(7),ij(13)},var12,var11),rr1,var11,I(1,2)),rr2,I(2,1),var21);
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,I(2,1),var21);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(7),ij(13)},var12,var11),rr1,I(1,1),I(1,2)),rr2,I(2,1),var21);
    end
    if includeL(14) && ~options.sep(6)
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,14},Zyb_st,Z2ba_sstt,gtr,ny,n2,nZyb,nZ2ba,rr2),rr2,var22,var21),var12,I(1,1),var11);
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,var22,var21);
    elseif includeL(14) && options.sep(6)
        %
        %
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,I(2,1),var21);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(7),ij(14)},var12,var11),rr1,I(1,1),var11),rr2,I(2,1),var21);
    end
    
    if includeL(12)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(constructPy2c(Q{7,12},Zyb_st,Z2ob_sst,gsr,ny,n2,nZyb,nZ2ob,rr2),rr2,I(2,1),var21);
        %Pop.Ry2{2} = Pop.Ry2{2} + int(constructPy2c(Q{7,12},Zyb_st,Z2ob_sst,gsr,ny,n2,nZyb,nZ2ob,rr2),rr2,I(2,1),var22);
        Pop.Ry2{2} = Pop.Ry2{2} + int(gsr * subs(N{ij(7),ij(12)},var12,var11),rr2,I(2,1),var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(gsr * subs(N{ij(7),ij(12)},var12,var11),rr2,I(2,1),var21);
    end
    if includeL(15) && ~options.sep(5)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{7,15},Zyb_st,Z2ab_sstt,gtr,ny,n2,nZyb,nZ2ab,rr2),rr2,I(2,1),var21),var12,var11,I(1,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,15},Zyb_st,Z2ab_sstt,gtr,ny,n2,nZyb,nZ2ab,rr2),rr2,I(2,1),var22),var12,var11,I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(15)},var12,var11),rr1,var11,I(1,2)),rr2,I(2,1),var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(7),ij(15)},var12,var11),rr1,var11,I(1,2)),rr2,I(2,1),var21);
    elseif includeL(15) && options.sep(5)==1
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{7,15},Zyb_st,Z2ab_sstt,gtr,ny,n2,nZyb,nZ2ab,rr2),rr2,I(2,1),var21),var12,I(1,1),I(1,2));
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,15},Zyb_st,Z2ab_sstt,gtr,ny,n2,nZyb,nZ2ab,rr2),rr2,I(2,1),var22),var12,I(1,1),I(1,2));
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(15)},var12,var11),rr1,I(1,1),I(1,2)),rr2,I(2,1),var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(7),ij(15)},var12,var11),rr1,I(1,1),I(1,2)),rr2,I(2,1),var21);
    end
    if includeL(16)
        %Pop.Ry2{3} = Pop.Ry2{3} + int(int(constructPy2c(Q{7,16},Zyb_st,Z2bb_sstt,gtr,ny,n2,nZyb,nZ2bb,rr2),rr2,I(2,1),var21),var12,I(1,1),var11);
        %Pop.Ry2{2} = Pop.Ry2{2} + int(int(constructPy2c(Q{7,16},Zyb_st,Z2bb_sstt,gtr,ny,n2,nZyb,nZ2bb,rr2),rr2,I(2,1),var22),var12,I(1,1),var11);
        Pop.Ry2{2} = Pop.Ry2{2} + int(int(grr * subs(N{ij(7),ij(16)},var12,var11),rr1,I(1,1),var11),rr2,I(2,1),var22);
        Pop.Ry2{3} = Pop.Ry2{3} + int(int(grr * subs(N{ij(7),ij(16)},var12,var11),rr1,I(1,1),var11),rr2,I(2,1),var21);
    end        
end

if options.sep(2)==1
    Pop.Ryy{3} = Pop.Ryy{2};
else
    Pop.Ryy{3} = var_swap(Pop.Ryy{2},var21,var22).';
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % Build R22
if includeL(8)
    %Pop.R22{1,1} = Pop.R22{1,1} + constructP22aa(Q{8,8},Z2oo_ss,Z2oo_ss,gss,n2,nZ2oo,nZ2oo);
    Pop.R22{1,1} = Pop.R22{1,1} + gss * subs(N{ij(8),ij(8)},var2,var1);
    
    if includeL(9) && ~options.sep(3)
        %Pop.R22{2,1} = Pop.R22{2,1} + constructP22aa(Q{8,9},Z2oo_ss,Z2ao_sst,gss,n2,nZ2oo,nZ2ao);
        Pop.R22{2,1} = Pop.R22{2,1} + gss * subs(N{ij(8),ij(9)},[rr1,var22],var1);
    elseif includeL(9) && options.sep(3)
        %Pop.R22{2,1} = Pop.R22{2,1} + constructP22aa(Q{8,9},Z2oo_ss,Z2ao_sst,gss,n2,nZ2oo,nZ2ao) + constructP22bb(Q{9,8},Z2ao_ss,Z2oo_sts,gts,n2,nZ2ao,nZ2oo,[1;0]);
        Pop.R22{2,1} = Pop.R22{2,1} + gss * subs(N{ij(8),ij(9)},[rr1,var22],var1) + gts * subs(N{ij(9),ij(8)},[rr1,var22],[var12;var21]);
    end
    if includeL(10)
        %Pop.R22{2,1} = Pop.R22{2,1} + constructP22bb(Q{10,8},Z2bo_sst,Z2oo_ss,gts,n2,nZ2bo,nZ2oo,[1;0]);
        Pop.R22{2,1} = Pop.R22{2,1} + gts * subs(N{ij(10),ij(8)},[rr1,var22],[var12;var21]);
    end
    
    if includeL(11) && ~options.sep(4)
        %Pop.R22{1,2} = Pop.R22{1,2} + constructP22aa(Q{8,11},Z2oo_ss,Z2oa_sst,gss,n2,nZ2oo,nZ2oa);
        Pop.R22{1,2} = Pop.R22{1,2} + gss * subs(N{ij(8),ij(11)},[var12,rr2],var1);
    elseif includeL(11) && options.sep(4)
        %Pop.R22{1,2} = Pop.R22{1,2} + constructP22aa(Q{8,11},Z2oo_ss,Z2oa_sst,gss,n2,nZ2oo,nZ2oa) + constructP22bb(Q{11,8},Z2oa_sst,Z2oo_ss,gst,n2,nZ2oa,nZ2oo,[0;1]);
        Pop.R22{1,2} = Pop.R22{1,2} + gss * subs(N{ij(8),ij(11)},[var12,rr2],var1) + gst * subs(N{ij(11),ij(8)},[var12,rr2],[var11;var22]);
    end
    if includeL(12)
        %Pop.R22{1,2} = Pop.R22{1,2} + constructP22bb(Q{12,8},Z2ob_sst,Z2oo_ss,gst,n2,nZ2ob,nZ2oo,[0;1]);
        Pop.R22{1,2} = Pop.R22{1,2} + gst * subs(N{ij(12),ij(8)},[var12,rr2],[var11;var22]);
    end
    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22aa(Q{8,13},Z2oo_ss,Z2aa_sstt,gss,n2,nZ2oo,nZ2aa);
        Pop.R22{2,2} = Pop.R22{2,2} + gss * subs(N{ij(8),ij(13)},rr,var1);
    elseif includeL(13) && ~options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + gss * subs(N{ij(8),ij(13)},rr,var1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gss * subs(N{ij(8),ij(13)},rr,var1);
    elseif includeL(13) && ~options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + gss * subs(N{ij(8),ij(13)},rr,var1);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gtt * subs(N{ij(13),ij(8)},rr,var2);
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22aa(Q{8,13},Z2oo_ss,Z2aa_sstt,gss,n2,nZ2oo,nZ2aa) +...
        %    constructP22bb(Q{13,8},Z2aa_sstt,Z2oo_ss,gtt,n2,nZ2aa,nZ2oo,[1;1]);
        Pop.R22{2,2} = Pop.R22{2,2} + gss * subs(N{ij(8),ij(13)},rr,var1) + gtt * subs(N{ij(13),ij(8)},rr,var2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gss * subs(N{ij(8),ij(13)},rr,var1) + gtt * subs(N{ij(13),ij(8)},rr,var2);
    end
    if includeL(14) && ~options.sep(6)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22aa(Q{8,14},Z2oo_ss,Z2ba_sstt,gss,n2,nZ2oo,nZ2ba);
        Pop.R22{3,2} = Pop.R22{3,2} + gss * subs(N{ij(8),ij(14)},rr,var1);
    elseif includeL(14) && options.sep(6)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22aa(Q{8,14},Z2oo_ss,Z2ba_stsst,gss,n2,nZ2oo,nZ2ba) +...
        %    constructP22bb(Q{14,8},Z2ba_stst,Z2oo_ss,gtt,n2,nZ2ba,nZ2oo,[1;1]);
        Pop.R22{2,2} = Pop.R22{2,2} + gtt * subs(N{ij(14),ij(8)},rr,var2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gss * subs(N{ij(8),ij(14)},rr,var1);
    end
    if includeL(15) && ~options.sep(5)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{15,8},Z2ab_sstt,Z2oo_ss,gtt,n2,nZ2ab,nZ2oo,[1;1]);
        Pop.R22{3,2} = Pop.R22{3,2} + gtt * subs(N{ij(15),ij(8)},rr,var2);
    elseif includeL(15) && options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + gtt * subs(N{ij(15),ij(8)},rr,var2);
        % % Pop.R22{3,2} = Pop.R22{3,2} + gtt * subs(N{ij(15),ij(8)},rr,var2);
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{16,8},Z2bb_sstt,Z2oo_ss,gtt,n2,nZ2bb,nZ2oo,[1;1]);
        Pop.R22{2,2} = Pop.R22{2,2} + gtt * subs(N{ij(16),ij(8)},rr,var2);
    end
    
end

if includeL(9) && ~options.sep(3)
    %Pop.R22{2,1} = Pop.R22{2,1} + int(constructP22cc(Q{9,9},Z2ao_sst,Z2ao_sst,grs,n2,nZ2ao,nZ2ao,[1;0],[rr1;rr2]),rr1,var11,I(1,2));
    Pop.R22{2,1} = Pop.R22{2,1} + int(grs * subs(N{ij(9),ij(9)},var22,var21),rr1,var11,I(1,2));
    if includeL(10)
        %Pop.R22{2,1} = Pop.R22{2,1} + int(constructP22cc(Q{10,9},Z2bo_sst,Z2ao_sst,grs,n2,nZ2bo,nZ2ao,[1;0],[rr1;rr2]),rr1,var12,var11);
        Pop.R22{2,1} = Pop.R22{2,1} + int(grs * subs(N{ij(10),ij(9)},var22,var21),rr1,var12,var11);
    end
    
    if includeL(11) && ~options.sep(4)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{9,11},Z2ao_sst,Z2oa_sst,gts,n2,nZ2ao,nZ2oa,[1;0]);
        Pop.R22{3,2} = Pop.R22{3,2} + gts * subs(N{ij(9),ij(11)},rr,[var12;var21]);
    elseif includeL(11) && options.sep(4)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{9,11},Z2ao_sst,Z2oa_sst,gts,n2,nZ2ao,nZ2oa,[1;0]);
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{11,9},Z2oa_sst,Z2ao_sst,gst,n2,nZ2oa,nZ2ao,[0;1]);
        Pop.R22{2,2} = Pop.R22{2,2} + gst * subs(N{ij(11),ij(9)},rr,[var11;var22]);
        Pop.R22{3,2} = Pop.R22{3,2} + gts * subs(N{ij(9),ij(11)},rr,[var12;var21]);
    end
    if includeL(12)
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{12,9},Z2ob_sst,Z2ao_sst,gst,n2,nZ2ob,nZ2ao,[0;1]);
        Pop.R22{2,2} = Pop.R22{2,2} + gst * subs(N{ij(12),ij(9)},rr,[var11;var22]);
    end
    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{9,13},Z2ao_sst,Z2aa_sstt,grs,n2,nZ2ao,nZ2aa,[1;0],[rr1;rr2]),rr1,var11,I(1,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{9,13},Z2ao_sst,Z2aa_sstt,grs,n2,nZ2ao,nZ2aa,[1;0],[rr1;rr2]),rr1,var12,I(1,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var11,I(1,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var12,I(1,2));
    elseif includeL(13) && ~options.sep(6)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var11,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var11,I(1,2));
    elseif includeL(13) && ~options.sep(5)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var11,I(1,2))...
                                    + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var11,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var12,I(1,2))...
        % %                             + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var12,I(1,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{9,13},Z2ao_sst,Z2aa_sstt,grs,n2,nZ2ao,nZ2aa,[1;0],[rr1;rr2]),rr1,var11,I(1,2)) +...
        %            int(constructP22cc(Q{13,9},Z2aa_sstt,Z2ao_sst,grt,n2,nZ2aa,nZ2ao,[1;1],[rr1;rr2]),rr1,var12,I(1,2));      
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var11,I(1,2))...
                                    + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var12,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var11,I(1,2))...
        % %                             + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var12,I(1,2));
    end
    if includeL(14) && ~options.sep(6)
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{9,14},Z2ao_sst,Z2ba_sstt,grs,n2,nZ2ao,nZ2ba,[1;0],[rr1;rr2]),rr1,var11,var12);
        Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(14)},rr2,var21),rr1,var11,var12);
    elseif includeL(14) && options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(14),ij(9)},rr2,var22),rr1,var12,var11);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(14)},rr2,var21),rr1,var11,var12);
    end
    if includeL(15) && ~options.sep(5)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{15,9},Z2ab_sstt,Z2ao_sst,grt,n2,nZ2ab,nZ2ao,[1;1],[rr1;rr2]),rr1,var11,I(1,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{15,9},Z2ab_sstt,Z2ao_sst,grt,n2,nZ2ab,nZ2ao,[1;1],[rr1;rr2]),rr1,var12,I(1,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,var11,I(1,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,var12,I(1,2));
    elseif includeL(15) && options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,var12,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,var12,I(1,2));
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{16,9},Z2bb_sstt,Z2ao_sst,grt,n2,nZ2bb,nZ2ao,[1;1],[rr1;rr2]),rr1,var12,var11);
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(16),ij(9)},rr2,var22),rr1,var12,var11);
    end
    
elseif includeL(9) && options.sep(3)
    %Pop.R22{2,1} = Pop.R22{2,1} + int(constructP22cc(Q{9,9},Z2ao_sst,Z2ao_sst,grs,n2,nZ2ao,nZ2ao,[1;0],[rr1;rr2]),rr1,I(1,1),I(1,2));
    Pop.R22{2,1} = Pop.R22{2,1} + int(grs * subs(N{ij(9),ij(9)},var22,var21),rr1,I(1,1),I(1,2));
    
    if includeL(11) && ~options.sep(4)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{9,11},Z2ao_stst,Z2oa_ss,gts,n2,nZ2ao,nZ2oa,[1;0]);
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{9,11},Z2ao_stst,Z2oa_ss,gts,n2,nZ2ao,nZ2oa,[1;0]);
        Pop.R22{2,2} = Pop.R22{2,2} + gts * subs(N{ij(9),ij(11)},rr,[var12;var21]);
        Pop.R22{3,2} = Pop.R22{3,2} + gts * subs(N{ij(9),ij(11)},rr,[var12;var21]);
    elseif includeL(11) && options.sep(4)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{9,11},Z2ao_sst,Z2oa_sst,gts,n2,nZ2ao,nZ2oa,[1;0]);
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{9,11},Z2ao_sst,Z2oa_sst,gts,n2,nZ2ao,nZ2oa,[1;0]);
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{11,9},Z2oa_sst,Z2ao_sst,gst,n2,nZ2oa,nZ2ao,[0;1]);
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{11,9},Z2oa_sst,Z2ao_sst,gst,n2,nZ2oa,nZ2ao,[0;1]);
        Pop.R22{2,2} = Pop.R22{2,2} + gst * subs(N{ij(11),ij(9)},rr,[var11;var22]) + gts * subs(N{ij(9),ij(11)},rr,[var12;var21]);
        Pop.R22{3,2} = Pop.R22{3,2} + gst * subs(N{ij(11),ij(9)},rr,[var11;var22]) + gts * subs(N{ij(9),ij(11)},rr,[var12;var21]);
    end
    if includeL(12)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{11,9},Z2ob_sst,Z2ao_sst,gst,n2,nZ2ob,nZ2ao,[0;1]);
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{12,9},Z2ob_sst,Z2ao_sst,gst,n2,nZ2ob,nZ2ao,[0;1]);
        Pop.R22{2,2} = Pop.R22{2,2} + gst * subs(N{ij(12),ij(9)},rr,[var11;var22]);
        Pop.R22{3,2} = Pop.R22{3,2} + gst * subs(N{ij(12),ij(9)},rr,[var11;var22]);
    end
    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{9,13},Z2ao_sst,Z2aa_sstt,grs,n2,nZ2ao,nZ2aa,[1;0],[rr1;rr2]),rr1,var12,I(1,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{9,13},Z2ao_sst,Z2aa_sstt,grs,n2,nZ2ao,nZ2aa,[1;0],[rr1;rr2]),rr1,var12,I(1,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var12,I(1,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var12,I(1,2));
    elseif includeL(13) && ~options.sep(6)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,I(1,1),I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,I(1,1),I(1,2));
    elseif includeL(13) && ~options.sep(5)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var12,I(1,2))...
                                    + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var11,I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,var12,I(1,2))...
        % %                             + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,var11,I(1,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{9,13},Z2ao_sst,Z2aa_sstt,grs,n2,nZ2ao,nZ2aa,[1;0],[rr1;rr2]),rr1,I(1,1),I(1,2)) +...
        %            int(constructP22cc(Q{13,9},Z2aa_sstt,Z2ao_sst,grt,n2,nZ2aa,nZ2ao,[1;1],[rr1;rr2]),rr1,I(1,1),I(1,2));      
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,I(1,1),I(1,2))...
                                    + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,I(1,1),I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(13)},rr2,var21),rr1,I(1,1),I(1,2))...
        % %                             + int(grt * subs(N{ij(13),ij(9)},rr2,var22),rr1,I(1,1),I(1,2));
    end
    if includeL(14) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{9,14},Z2ao_sst,Z2ba_sstt,grs,n2,nZ2ao,nZ2ba,[1;0],[rr1;rr2]),rr1,I(1,1),var12);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{9,14},Z2ao_sst,Z2ba_sstt,grs,n2,nZ2ao,nZ2ba,[1;0],[rr1;rr2]),rr1,I(1,1),var12);
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(14)},rr2,var21),rr1,I(1,1),var12);
        Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(14)},rr2,var21),rr1,I(1,1),var12);
    elseif includeL(14) && options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(9),ij(14)},rr2,var21),rr1,I(1,1),var12)...
                                    + int(grt * subs(N{ij(14),ij(9)},rr2,var22),rr1,I(1,1),var11);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(9),ij(14)},rr2,var21),rr1,I(1,1),var12)...
        % %                             + int(grt * subs(N{ij(14),ij(9)},rr2,var22),rr1,I(1,1),var11);
    end
    if includeL(15) && ~options.sep(5)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{15,9},Z2ab_sstt,Z2ao_sst,grt,n2,nZ2ab,nZ2ao,[1;1],[rr1;rr2]),rr1,var11,I(1,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{15,9},Z2ab_sstt,Z2ao_sst,grt,n2,nZ2ab,nZ2ao,[1;1],[rr1;rr2]),rr1,var11,I(1,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,var11,I(1,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,var11,I(1,2));
    elseif includeL(15) && options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,I(1,1),I(1,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(9)},rr2,var22),rr1,I(1,1),I(1,2));
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{16,9},Z2bb_sstt,Z2ao_sst,grt,n2,nZ2bb,nZ2ao,[1;1],[rr1;rr2]),rr1,I(1,1),var11);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{16,9},Z2bb_sstt,Z2ao_sst,grt,n2,nZ2bb,nZ2ao,[1;1],[rr1;rr2]),rr1,I(1,1),var11);
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(16),ij(9)},rr2,var22),rr1,I(1,1),var11);
        Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(16),ij(9)},rr2,var22),rr1,I(1,1),var11);
    end
    
end

if includeL(10)
    %Pop.R22{2,1} = Pop.R22{2,1} + int(constructP22cc(Q{10,10},Z2bo_sst,Z2bo_sst,grs,n2,nZ2bo,nZ2bo,[1;0],[rr1;rr2]),rr1,I(1,1),var12);
    Pop.R22{2,1} = Pop.R22{2,1} + int(grs * subs(N{ij(10),ij(10)},var22,var21),rr1,I(1,1),var12);
    
    if includeL(11) && ~options.sep(4)
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{10,11},Z2bo_sst,Z2oa_sst,gts,n2,nZ2bo,nZ2oa,[1;0]);
        Pop.R22{2,2} = Pop.R22{2,2} + gts * subs(N{ij(10),ij(11)},rr,[var12;var21]);
    elseif includeL(11) && options.sep(4)
        %Pop.R22{2,2} = Pop.R22{2,2} + constructP22bb(Q{10,11},Z2bo_sst,Z2oa_sst,gts,n2,nZ2bo,nZ2oa,[1;0]);
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{11,10},Z2ob_sst,Z2bo_sst,gst,n2,nZ2ob,nZ2bo,[0;1]);
        Pop.R22{2,2} = Pop.R22{2,2} + gts * subs(N{ij(10),ij(11)},rr,[var12;var21]);
        Pop.R22{3,2} = Pop.R22{3,2} + gst * subs(N{ij(11),ij(10)},rr,[var11;var22]);
    end
    if includeL(12)
        %Pop.R22{3,2} = Pop.R22{3,2} + constructP22bb(Q{12,10},Z2ob_sst,Z2bo_sst,gst,n2,nZ2ob,nZ2bo,[0;1]);
        Pop.R22{3,2} = Pop.R22{3,2} + gst * subs(N{ij(12),ij(10)},rr,[var11;var22]);
    end
    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{10,13},Z2bo_sst,Z2aa_sstt,grs,n2,nZ2bo,nZ2aa,[1;0],[rr1;rr2]),rr1,var12,var11);
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(10),ij(13)},rr2,var21),rr1,var12,var11);
    elseif includeL(13) && ~options.sep(6)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(10),ij(13)},rr2,var21),rr1,I(1,1),var11);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(10),ij(13)},rr2,var21),rr1,I(1,1),var11);
    elseif includeL(13) && ~options.sep(5)
        %
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(10),ij(13)},rr2,var21),rr1,var12,var11);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(13),ij(10)},rr2,var22),rr1,var11,var12);
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{10,13},Z2bo_sst,Z2aa_sstt,grs,n2,nZ2bo,nZ2aa,[1;0],[rr1;rr2]),rr1,I(1,1),var11) +...
        %            int(constructP22cc(Q{13,9},Z2aa_sstt,Z2ao_sst,grt,n2,nZ2aa,nZ2ao,[1;1],[rr1;rr2]),rr1,I(1,1),var12);      
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(10),ij(13)},rr2,var21),rr1,I(1,1),var11)...
                                    + int(grt * subs(N{ij(13),ij(10)},rr2,var22),rr1,I(1,1),var12);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(10),ij(13)},rr2,var21),rr1,I(1,1),var11)...
        % %                             + int(grt * subs(N{ij(13),ij(10)},rr2,var22),rr1,I(1,1),var12);
    end
    if includeL(14) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{10,14},Z2bo_sst,Z2ba_sstt,grs,n2,nZ2bo,nZ2ba,[1;0],[rr1;rr2]),rr1,I(1,1),var12);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{10,14},Z2bo_sst,Z2ba_sstt,grs,n2,nZ2bo,nZ2ba,[1;0],[rr1;rr2]),rr1,I(1,1),var11);
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(10),ij(14)},rr2,var21),rr1,I(1,1),var12);
        Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(10),ij(14)},rr2,var21),rr1,I(1,1),var11);
    elseif includeL(14) && options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grs * subs(N{ij(10),ij(14)},rr2,var21),rr1,I(1,1),var12)...
                                    + int(grt * subs(N{ij(14),ij(10)},rr2,var22),rr1,I(1,1),var12);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grs * subs(N{ij(10),ij(14)},rr2,var21),rr1,I(1,1),var11)...
        % %                             + int(grt * subs(N{ij(14),ij(10)},rr2,var22),rr1,I(1,1),var11);
    end
    if includeL(15) && ~options.sep(5)
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{15,10},Z2ab_sstt,Z2bo_sst,grt,n2,nZ2ab,nZ2bo,[1;1],[rr1;rr2]),rr1,var11,var12);
        Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(10)},rr2,var22),rr1,var11,var12);
    elseif includeL(15) && options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(15),ij(10)},rr2,var22),rr1,I(1,1),var12);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(15),ij(10)},rr2,var22),rr1,I(1,1),var12);
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{16,10},Z2bb_sstt,Z2bo_sst,grt,n2,nZ2bb,nZ2bo,[1;1],[rr1;rr2]),rr1,I(1,1),var12);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{16,10},Z2bb_sstt,Z2bo_sst,grt,n2,nZ2bb,nZ2bo,[1;1],[rr1;rr2]),rr1,I(1,1),var11);
        Pop.R22{2,2} = Pop.R22{2,2} + int(grt * subs(N{ij(16),ij(10)},rr2,var22),rr1,I(1,1),var12);
        Pop.R22{3,2} = Pop.R22{3,2} + int(grt * subs(N{ij(16),ij(10)},rr2,var22),rr1,I(1,1),var11);
    end
    
end

if includeL(11) && ~options.sep(4)
    %Pop.R22{1,2} = Pop.R22{1,2} + int(constructP22cc(Q{11,11},Z2oa_sst,Z2oa_sst,gsr,n2,nZ2oa,nZ2oa,[0;1],[rr1;rr2]),rr2,var21,I(2,2));
    Pop.R22{1,2} = Pop.R22{1,2} + int(gsr * subs(N{ij(11),ij(11)},var12,var11),rr2,var21,I(2,2));
    if includeL(12)
        %Pop.R22{1,2} = Pop.R22{1,2} + int(constructP22cc(Q{12,11},Z2ob_sst,Z2oa_sst,gsr,n2,nZ2ob,nZ2oa,[0;1],[rr1;rr2]),rr2,var22,var21);
        Pop.R22{1,2} = Pop.R22{1,2} + int(gsr * subs(N{ij(12),ij(11)},var12,var11),rr2,var22,var21);
    end
    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{11,13},Z2oa_sst,Z2aa_sstt,gsr,n2,nZ2oa,nZ2aa,[0;1],[rr1;rr2]),rr2,var21,I(2,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{13,11},Z2aa_sstt,Z2oa_sst,gtr,n2,nZ2aa,nZ2oa,[1;1],[rr1;rr2]),rr2,var21,I(2,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var21,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var21,I(2,2))...
                                    + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var21,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var21,I(2,2))...
        % %                             + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var21,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var22,I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{11,13},Z2oa_sst,Z2aa_sstt,gsr,n2,nZ2oa,nZ2aa,[0;1],[rr1;rr2]),rr2,var21,I(2,2)) +...
        %            int(constructP22cc(Q{13,11},Z2aa_sstt,Z2oa_sst,gtr,n2,nZ2aa,nZ2oa,[1;1],[rr1;rr2]),rr2,var22,I(2,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var21,I(2,2))...
                                    + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var22,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var22,I(2,2))...
        % %                             + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var21,I(2,2));
    end
    if includeL(14) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{14,11},Z2ba_sstt,Z2oa_sst,gtr,n2,nZ2ba,nZ2oa,[1;1],[rr1;rr2]),rr2,var21,I(2,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{11,14},Z2oa_sst,Z2ba_sstt,gsr,n2,nZ2oa,nZ2ba,[0;1],[rr1;rr2]),rr2,var21,I(2,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(14),ij(11)},rr1,var12),rr2,var21,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(11),ij(14)},rr1,var11),rr2,var21,I(2,2));
    elseif includeL(14) && options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(14),ij(11)},rr1,var12),rr2,var22,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(11),ij(14)},rr1,var11),rr2,var21,I(2,2));
    end
    if includeL(15) && ~options.sep(5)
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{15,11},Z2ab_sstt,Z2oa_sst,gtr,n2,nZ2ab,nZ2oa,[1;1],[rr1;rr2]),rr2,var22,var21);
        Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(11)},rr1,var12),rr2,var22,var21);
    elseif includeL(15) && options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(15),ij(11)},rr1,var12),rr2,var22,var21);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(11)},rr1,var12),rr2,var22,var21);
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{16,11},Z2bb_sstt,Z2oa_sst,gtr,n2,nZ2bb,nZ2oa,[1;1],[rr1;rr2]),rr2,var22,var21);
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(16),ij(11)},rr1,var12),rr2,var22,var21);
    end
    
elseif includeL(11) && options.sep(4)
    %Pop.R22{1,2} = Pop.R22{1,2} + int(constructP22cc(Q{11,11},Z2oa_sst,Z2oa_sst,gsr,n2,nZ2oa,nZ2oa,[0;1],[rr1;rr2]),rr2,I(2,1),I(2,2));
    Pop.R22{1,2} = Pop.R22{1,2} + int(gsr * subs(N{ij(11),ij(11)},var12,var11),rr2,I(2,1),I(2,2));
    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{11,13},Z2oa_sst,Z2aa_sstt,gsr,n2,nZ2oa,nZ2aa,[0;1],[rr1;rr2]),rr2,var22,I(2,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{13,11},Z2aa_sst,Z2oa_stst,gtr,n2,nZ2aa,nZ2oa,[1;1],[rr1;rr2]),rr2,var21,I(2,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var22,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var21,I(2,2));
    elseif includeL(13) && ~options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var22,I(2,2))...
                                    + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var21,I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,var21,I(2,2))...
        % %                             + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,var22,I(2,2));
    elseif includeL(13) && ~options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,I(2,1),I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,I(2,1),I(2,2));
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{11,13},Z2oa_sst,Z2aa_sstt,gsr,n2,nZ2oa,nZ2aa,[0;1],[rr1;rr2]),rr2,I(2,1),I(2,2)) +...
        %            int(constructP22cc(Q{13,11},Z2aa_sstt,Z2oa_sst,gtr,n2,nZ2aa,nZ2oa,[1;1],[rr1;rr2]),rr2,I(2,1),I(2,2));      
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,I(2,1),I(2,2))...
                                    + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,I(2,1),I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(11)},rr1,var12),rr2,I(2,1),I(2,2))...
        % %                             + int(gsr * subs(N{ij(11),ij(13)},rr1,var11),rr2,I(2,1),I(2,2));
    end    
    if includeL(14) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{14,11},Z2ba_sstt,Z2oa_sst,gtr,n2,nZ2ba,nZ2oa,[1;1],[rr1;rr2]),rr2,var21,I(2,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{11,14},Z2oa_sst,Z2ba_sstt,gsr,n2,nZ2oa,nZ2ba,[0;1],[rr1;rr2]),rr2,var22,I(2,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(14),ij(11)},rr1,var12),rr2,var21,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(11),ij(14)},rr1,var11),rr2,var22,I(2,2));
    elseif includeL(14) && options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(14),ij(11)},rr1,var12),rr2,I(2,1),I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(11),ij(14)},rr1,var11),rr2,I(2,1),I(2,2));
    end
    if includeL(15) && ~options.sep(5)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{11,15},Z2oa_sst,Z2ab_sstt,gsr,n2,nZ2oa,nZ2ab,[0;1],[rr1;rr2]),rr2,I(2,1),var22);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{15,11},Z2ab_sst,Z2oa_stst,gtr,n2,nZ2ab,nZ2oa,[1;1],[rr1;rr2]),rr2,I(2,1),var21);
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(15)},rr1,var11),rr2,I(2,1),var22);
        Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(11)},rr1,var12),rr2,I(2,1),var21);
    elseif includeL(15) && options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(11),ij(15)},rr1,var11),rr2,I(2,1),var22)...
                                    + int(gtr * subs(N{ij(15),ij(11)},rr1,var12),rr2,I(2,1),var21);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(11)},rr1,var12),rr2,I(2,1),var21)...
        % %                             + int(gsr * subs(N{ij(11),ij(15)},rr1,var11),rr2,I(2,1),var22);
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{16,11},Z2bb_sstt,Z2oa_sst,gtr,n2,nZ2bb,nZ2oa,[1;1],[rr1;rr2]),rr2,I(2,1),var21);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{16,11},Z2bb_sstt,Z2oa_sst,gtr,n2,nZ2bb,nZ2oa,[1;1],[rr1;rr2]),rr2,I(2,1),var22);
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(16),ij(11)},rr1,var12),rr2,I(2,1),var21);
        Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(11),ij(16)},rr1,var11),rr2,I(2,1),var22);
    end
    
end

if includeL(12)
    %Pop.R22{1,2} = Pop.R22{1,2} + int(constructP22cc(Q{12,12},Z2ob_sst,Z2ob_sst,gsr,n2,nZ2ob,nZ2ob,[0;1],[rr1;rr2]),rr2,I(2,1),var22);
    Pop.R22{1,2} = Pop.R22{1,2} + int(gsr * subs(N{ij(12),ij(12)},var12,var11),rr2,I(2,1),var22);
    
    if includeL(13) && ~options.sep(5) && ~options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{12,13},Z2ob_sst,Z2aa_sstt,gsr,n2,nZ2oa,nZ2aa,[0;1],[rr1;rr2]),rr2,var22,var21);
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(12),ij(13)},rr1,var11),rr2,var22,var21);
    elseif includeL(13) && ~options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(12),ij(13)},rr1,var11),rr2,var22,var21);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(12),ij(13)},rr1,var11),rr2,var22,var21);
    elseif includeL(13) && ~options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(12),ij(13)},rr1,var11),rr2,I(2,1),var21);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(12)},rr1,var12),rr2,I(2,1),var22);
    elseif includeL(13) && options.sep(5) && options.sep(6)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{12,13},Z2ob_sst,Z2aa_sstt,gsr,n2,nZ2ob,nZ2aa,[0;1],[rr1;rr2]),rr2,I(2,1),var21) +...
        %            int(constructP22cc(Q{13,12},Z2aa_sstt,Z2ob_sst,gtr,n2,nZ2aa,nZ2ob,[1;1],[rr1;rr2]),rr2,I(2,1),var22);
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(12),ij(13)},rr1,var11),rr2,I(2,1),var21)...
                                    + int(gtr * subs(N{ij(13),ij(12)},rr1,var12),rr2,I(2,1),var22);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(13),ij(12)},rr1,var12),rr2,I(2,1),var22)...
        % %                             + int(gsr * subs(N{ij(12),ij(13)},rr1,var11),rr2,I(2,1),var21);
    end 
    
    if includeL(14) && ~options.sep(6)
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{12,14},Z2ob_sst,Z2ba_sstt,gsr,n2,nZ2ob,nZ2ba,[0;1],[rr1;rr2]),rr2,var22,var21);
        Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(12),ij(14)},rr1,var11),rr2,var22,var21);
    elseif includeL(14) && options.sep(6)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(14),ij(12)},rr1,var12),rr2,I(2,1),var22);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(12),ij(14)},rr1,var11),rr2,I(2,1),var21);
    end
    if includeL(15) && ~options.sep(5)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{12,15},Z2ob_sst,Z2ab_sstt,gsr,n2,nZ2ob,nZ2ab,[0;1],[rr1;rr2]),rr2,I(2,1),var22);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{15,12},Z2ab_sstt,Z2ob_sst,gtr,n2,nZ2ab,nZ2ob,[1;1],[rr1;rr2]),rr2,I(2,1),var22);
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(12),ij(15)},rr1,var11),rr2,I(2,1),var22);
        Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(12)},rr1,var12),rr2,I(2,1),var22);
    elseif includeL(15) && options.sep(5)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(gsr * subs(N{ij(12),ij(15)},rr1,var11),rr2,I(2,1),var22)...
                                    + int(gtr * subs(N{ij(15),ij(12)},rr1,var12),rr2,I(2,1),var22);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(gtr * subs(N{ij(15),ij(12)},rr1,var12),rr2,I(2,1),var22)...
        % %                             + int(gsr * subs(N{ij(12),ij(15)},rr1,var11),rr2,I(2,1),var22);
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(constructP22cc(Q{16,12},Z2bb_sstt,Z2ob_sst,gtr,n2,nZ2bb,nZ2ob,[1;1],[rr1;rr2]),rr2,I(2,1),var22);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(constructP22cc(Q{12,16},Z2ob_sst,Z2bb_sstt,gsr,n2,nZ2ob,nZ2bb,[0;1],[rr1;rr2]),rr2,I(2,1),var22);
        Pop.R22{2,2} = Pop.R22{2,2} + int(gtr * subs(N{ij(16),ij(12)},rr1,var12),rr2,I(2,1),var22);
        Pop.R22{3,2} = Pop.R22{3,2} + int(gsr * subs(N{ij(12),ij(16)},rr1,var11),rr2,I(2,1),var22);
    end    
end

if includeL(13) && ~options.sep(5) && ~options.sep(6)
    %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{13,13},Z2aa_sstt,Z2aa_sstt,grr,n2,nZ2aa,nZ2aa,[1;1],[rr1;rr2]),rr1,var11,I(1,2)),rr2,var21,I(2,2));
    %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{13,13},Z2aa_sstt,Z2aa_sstt,grr,n2,nZ2aa,nZ2aa,[1;1],[rr1;rr2]),rr1,var12,I(1,2)),rr2,var21,I(2,2));
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(13),ij(13)},rr1,var11,I(1,2)),rr2,var21,I(2,2));
    Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(13)},rr1,var12,I(1,2)),rr2,var21,I(2,2));
    
    if includeL(14)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{14,13},Z2ba_sstt,Z2aa_sstt,grr,n2,nZ2ba,nZ2aa,[1;1],[rr1;rr2]),rr1,var12,var11),rr2,var21,I(2,2));
        %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{13,14},Z2aa_sstt,Z2ba_sstt,grr,n2,nZ2aa,nZ2ba,[1;1],[rr1;rr2]),rr1,var11,var12),rr2,var21,I(2,2));
        Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(14),ij(13)},rr1,var12,var11),rr2,var21,I(2,2));
        Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(14)},rr1,var11,var12),rr2,var21,I(2,2));
    end
    if includeL(15)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{15,13},Z2ab_sstt,Z2aa_sstt,grr,n2,nZ2ab,nZ2aa,[1;1],[rr1;rr2]),rr1,var11,I(1,2)),rr2,var22,var21);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{15,13},Z2ab_sstt,Z2aa_sstt,grr,n2,nZ2ab,nZ2aa,[1;1],[rr1;rr2]),rr1,var12,I(1,2)),rr2,var22,var21);
        Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(15),ij(13)},rr1,var11,I(1,2)),rr2,var22,var21);
        Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(13)},rr1,var12,I(1,2)),rr2,var22,var21);
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{16,13},Z2bb_sstt,Z2aa_sstt,grr,n2,nZ2bb,nZ2aa,[1;1],[rr1;rr2]),rr1,var12,var11),rr2,var22,var21);
        Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(16),ij(13)},rr1,var12,var11),rr2,var22,var21);
    end
    
elseif includeL(13) && ~options.sep(6)
    %
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(13),ij(13)},rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(13)},rr1,I(1,1),I(1,2)),rr2,var21,I(2,2));
    
    if includeL(15)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(15),ij(13)},rr1,I(1,1),I(1,2)),rr2,var22,var21);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(13)},rr1,I(1,1),I(1,2)),rr2,var22,var21);
    end
    
elseif includeL(13) && ~options.sep(5)
    %
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(13),ij(13)},rr1,var11,I(1,2)),rr2,I(2,1),I(2,2));
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(13)},rr1,var12,I(1,2)),rr2,I(2,1),I(2,2));
    
    if includeL(14)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(14),ij(13)},rr1,var12,var11),rr2,I(2,1),I(2,2));
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(14)},rr1,var11,var12),rr2,I(2,1),I(2,2));
    end

elseif includeL(13) && options.sep(5) && options.sep(6)
    %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{13,13},Z2aa_sstt,Z2aa_sstt,grr,n2,nZ2aa,nZ2aa,[1;1],[rr1;rr2]),rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(13),ij(13)},rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(13),ij(13)},rr1,I(1,1),I(1,2)),rr2,I(2,1),I(2,2));
end

if includeL(14) && ~options.sep(6)
    %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{14,14},Z2ba_sstt,Z2ba_sstt,grr,n2,nZ2ba,nZ2ba,[1;1],[rr1;rr2]),rr1,I(1,1),var12),rr2,var21,I(2,2));
    %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{14,14},Z2ba_sstt,Z2ba_sstt,grr,n2,nZ2ba,nZ2ba,[1;1],[rr1;rr2]),rr1,I(1,1),var11),rr2,var21,I(2,2));
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(14),ij(14)},rr1,I(1,1),var12),rr2,var21,I(2,2));
    Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(14),ij(14)},rr1,I(1,1),var11),rr2,var21,I(2,2));
    
    if includeL(15)
        %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{15,14},Z2ab_sstt,Z2ba_sstt,grr,n2,nZ2ab,nZ2ba,[1;1],[rr1;rr2]),rr1,var11,var12),rr2,var22,var21);
        Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(14)},rr1,var11,var12),rr2,var22,var21);
    end
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{16,14},Z2bb_sstt,Z2ba_sstt,grr,n2,nZ2bb,nZ2ba,[1;1],[rr1;rr2]),rr1,I(1,1),var12),rr2,var22,var21);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{16,14},Z2bb_sstt,Z2ba_sstt,grr,n2,nZ2bb,nZ2ba,[1;1],[rr1;rr2]),rr1,I(1,1),var11),rr2,var22,var21);
        Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(16),ij(14)},rr1,I(1,1),var12),rr2,var22,var21);
        Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(16),ij(14)},rr1,I(1,1),var11),rr2,var22,var21);
    end
    
elseif includeL(14) && options.sep(6)
    %
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(14),ij(14)},rr1,I(1,1),var12),rr2,I(2,1),I(2,2));
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(14),ij(14)},rr1,I(1,1),var11),rr2,I(2,1),I(2,2));
    
    if includeL(15)
        %
        Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(14),ij(15)},rr1,var12,var11),rr2,I(2,1),var22);
        % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(14)},rr1,var11,var12),rr2,I(2,1),var21);
    end
    
end

if includeL(15) && ~options.sep(5)
    %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{15,15},Z2ab_sstt,Z2ab_sstt,grr,n2,nZ2ab,nZ2ab,[1;1],[rr1;rr2]),rr1,var11,I(1,2)),rr2,I(2,1),var22);
    %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{15,15},Z2ab_sstt,Z2ab_sstt,grr,n2,nZ2ab,nZ2ab,[1;1],[rr1;rr2]),rr1,var12,I(1,2)),rr2,I(2,1),var22);
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(15),ij(15)},rr1,var11,I(1,2)),rr2,I(2,1),var22);
    Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(15)},rr1,var12,I(1,2)),rr2,I(2,1),var22);
    
    if includeL(16)
        %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{16,15},Z2bb_sstt,Z2ab_sstt,grr,n2,nZ2bb,nZ2ab,[1;1],[rr1;rr2]),rr1,var12,var11),rr2,I(2,1),var22);
        %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{15,16},Z2ab_sstt,Z2bb_sstt,grr,n2,nZ2ab,nZ2bb,[1;1],[rr1;rr2]),rr1,var11,var12),rr2,I(2,1),var22);
        Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(16),ij(15)},rr1,var12,var11),rr2,I(2,1),var22);
        Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(16)},rr1,var11,var12),rr2,I(2,1),var22);
    end
    
elseif includeL(15) && options.sep(5)
    %
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(15),ij(15)},rr1,I(1,1),I(1,2)),rr2,I(2,1),var22);
    % % Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(15),ij(15)},rr1,I(1,1),I(1,2)),rr2,I(2,1),var22);
    
end

if includeL(16)
    %Pop.R22{2,2} = Pop.R22{2,2} + int(int(constructP22cc(Q{16,16},Z2bb_sstt,Z2bb_sstt,grr,n2,nZ2bb,nZ2bb,[1;1],[rr1;rr2]),rr1,I(1,1),var12),rr2,I(2,1),var22);
    %Pop.R22{3,2} = Pop.R22{3,2} + int(int(constructP22cc(Q{16,16},Z2bb_sstt,Z2bb_sstt,grr,n2,nZ2bb,nZ2bb,[1;1],[rr1;rr2]),rr1,I(1,1),var11),rr2,I(2,1),var22);
    Pop.R22{2,2} = Pop.R22{2,2} + int(int(grr * N{ij(16),ij(16)},rr1,I(1,1),var12),rr2,I(2,1),var22);
    Pop.R22{3,2} = Pop.R22{3,2} + int(int(grr * N{ij(16),ij(16)},rr1,I(1,1),var11),rr2,I(2,1),var22);
end

if options.sep(3)==1
    Pop.R22{3,1} = Pop.R22{2,1};
else
    Pop.R22{3,1} = var_swap(Pop.R22{2,1},var11,var12).';
end
if options.sep(4)==1
    Pop.R22{1,3} = Pop.R22{1,2};
else
    Pop.R22{1,3} = var_swap(Pop.R22{1,2},var21,var22).';
end
if options.sep(5) && options.sep(6)
    Pop.R22{3,2} = Pop.R22{2,2};
    Pop.R22{2,3} = Pop.R22{2,2};
    Pop.R22{3,3} = Pop.R22{2,2};
elseif options.sep(5)
    Pop.R22{3,2} = Pop.R22{2,2};
    Pop.R22{2,3} = var_swap(var_swap(Pop.R22{3,2},var11,var12),var21,var22).';
    Pop.R22{3,3} = Pop.R22{2,3};
elseif options.sep(6)
    Pop.R22{2,3} = Pop.R22{2,2};
    Pop.R22{3,2} = var_swap(var_swap(Pop.R22{2,3},var11,var12),var21,var22).';
    Pop.R22{3,3} = Pop.R22{3,2};
else
    Pop.R22{3,3} = var_swap(var_swap(Pop.R22{2,2},var11,var12),var21,var22).';
    Pop.R22{2,3} = var_swap(var_swap(Pop.R22{3,2},var11,var12),var21,var22).';
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

Pop.Rx0 = Pop.R0x.';
Pop.Ry0 = Pop.R0y.';
Pop.R20 = Pop.R02.';

Pop.Ryx = Pop.Rxy.';
Pop.R2x{1,1} = Pop.Rx2{1,1}.';
Pop.R2x{2,1} = var_swap(Pop.Rx2{3,1}.',var11,var12);
Pop.R2x{3,1} = var_swap(Pop.Rx2{2,1}.',var11,var12);
Pop.R2y{1,1} = Pop.Ry2{1,1}.';
Pop.R2y{1,2} = var_swap(Pop.Ry2{1,3}.',var21,var22);
Pop.R2y{1,3} = var_swap(Pop.Ry2{1,2}.',var21,var22);


end
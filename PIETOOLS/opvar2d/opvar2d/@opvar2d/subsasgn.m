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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07_06_2021

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
        error('()- like subsassign is currently not supported for opvar2d objects.');
        
%         if length(L)==1
%             temp = polynomial(RHS);
%         else
%             % Peform all subsasgn but L(1)
%             temp = subsref(a,L(1));
%             temp = subsasgn(temp,L(2:end),RHS);
%         end
%         
%         %  Three '()'-subsasgn cases
%         if length(L(1).subs)==1 && strcmp(L(1).subs{1},':')
%             b = PVsubsasgn_colon(a,L,temp);
%         elseif length(L(1).subs)==1
%             b = PVsubsasgn_1idx(a,L,temp);
%         else
%             if strcmp(L(1).subs{1},':')
%                 % To handle the case a(:,1)=x1 when a is undefined or 0-by-0
%                 if all(sza==[0 0])
%                     sza(1) = 1;
%                 end
%                 L(1).subs{1} = 1:sza(1);
%             end
%             if strcmp(L(1).subs{2},':')
%                 % To handle the case a(1,:)=x1 when a is undefined or 0-by-0
%                 if all(sza==[0 0])
%                     sza(2) = 1;
%                 end
%                 L(1).subs{2} = 1:sza(2);
%             end
%             b = PVsubsasgn_2idx(a,L,temp);
%         end
        
    case '{}'
        error('{}- like subsassign is not supported for opvar2d objects.');
end

b.dim = b.dim;

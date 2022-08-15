function [Pplus] = plus(P1,P2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pplus] = plus(P1,P2) performs addition of two operators 
%   P1,P2: R^m0 x L2^mx x L2^my x L2^m2 to R^n0 x L2^nx x L2^ny x L2^n2
% Date: 01/31/21
% Version: 1.0
% 
% INPUT
% P1, P2: opvar2d class objects
% 
% OUTPUT
% Pplus: returns P1+P2
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - plus
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
% Initial coding DJ - 01_31_2021  
%   ^ Based heavily on "@opvar"-plus code by SS ^



if (isa(P1,'polynomial') || isa(P1,'double')) && all(size(P1)==[1,1]) && all(P2.dim(:,1)==P2.dim(:,2))
    Pplus = P2;
    Pplus.R00 = Pplus.R00 + P1*eye(P2.dim(1,2));
    Pplus.Rxx{1,1} = Pplus.Rxx{1,1} + P1*eye(P2.dim(2,2));
    Pplus.Ryy{1,1} = Pplus.Ryy{1,1} + P1*eye(P2.dim(3,2));
    Pplus.R22{1,1} = Pplus.R22{1,1} + P1*eye(P2.dim(4,2));
    
elseif (isa(P2,'polynomial') || isa(P2,'double')) && all(size(P2)==[1,1]) && all(P1.dim(:,1)==P1.dim(:,2))
    Pplus = P1;
    Pplus.R00 = Pplus.R00 + P2*eye(P1.dim(1,2));
    Pplus.Rxx{1,1} = Pplus.Rxx{1,1} + P2*eye(P1.dim(2,2));
    Pplus.Ryy{1,1} = Pplus.Ryy{1,1} + P2*eye(P1.dim(3,2));
    Pplus.R22{1,1} = Pplus.R22{1,1} + P2*eye(P1.dim(4,2));
    
elseif ~isa(P1,'opvar2d')
    P2dim = P2.dim;
    if all(all(P2dim==0)) || isempty(P1)
        Pplus = P1;
    elseif all(all(P2dim(2:4,:)==zeros(3,2))) && all(size(P1)==P2dim(1,:))
        Pplus = P2;
        Pplus.R00 = P1 + Pplus.R00;
    elseif all(all(P2dim([1,3:4],:)==zeros(3,2))) && all(size(P1)==P2dim(2,:))
        Pplus = P2;
        Pplus.Rxx{1,1} = P1 + Pplus.Rxx{1,1};
    elseif all(all(P2dim([1:2,4],:)==zeros(3,2))) && all(size(P1)==P2dim(3,:))
        Pplus = P2;
        Pplus.Ryy{1,1} = P1 + Pplus.Ryy{1,1};
    elseif all(all(P2dim(1:3,:)==zeros(3,2))) && all(size(P1)==P2dim(4,:))
        Pplus = P2;
        Pplus.R22{1,1} = P1 + Pplus.R22{1,1};
    end
    
elseif ~isa(P2,'opvar2d')
    P1dim = P1.dim;
    if all(all(P1dim==0)) || isempty(P2)
        Pplus = P2;
    elseif all(all(P1dim(2:4,:)==zeros(3,2))) && all(size(P2)==P1dim(1,:))
        Pplus = P1;
        Pplus.R00 = P2 + Pplus.R00;
    elseif all(all(P1dim([1,3:4],:)==zeros(3,2))) && all(size(P2)==P1dim(2,:))
        Pplus = P1;
        Pplus.Rxx{1,1} = P2 + Pplus.Rxx{1,1};
    elseif all(all(P1dim([1:2,4],:)==zeros(3,2))) && all(size(P2)==P1dim(3,:))
        Pplus = P1;
        Pplus.Ryy{1,1} = P2 + Pplus.Ryy{1,1};
    elseif all(all(P1dim(1:3,:)==zeros(3,2))) && all(size(P2)==P1dim(4,:))
        Pplus = P1;
        Pplus.R22{1,1} = P2 + Pplus.R22{1,1};
    end
    
else %both are PI operators
    if any(any(P1.dim~=P2.dim))
        error('Operator dimensions do not match. Check the validity of the addition operation.');
    end
    if any(P1.I~=P2.I)
        error('Operators act on different intervals and cannot be added');
    end
    %Create a holder variable for the resultant operator
    Pplus = P1;
    
    fset = {'R00', 'R0x', 'R0y', 'R02', 'Rx0', 'Rxy', 'Ry0', 'Ryx', 'R20'};
    for f=fset
        if ~isempty(P1.(f{:})) && ~isempty(P2.(f{:}))
            Pplus.(f{:}) = P1.(f{:}) + P2.(f{:});
        elseif isempty(P1.(f{:}))
            Pplus.(f{:}) = P2.(f{:});
        elseif isempty(P2.(f{:}))
            Pplus.(f{:}) = P1.(f{:});
        else
            Pplus.(f{:}) = [];
        end
    end 
    
    fset = {'Rxx','Rx2','R2x'};
    for f=fset
        P1R = P1.(f{:});    P2R = P2.(f{:});
        PplusR = cell(3,1);
        for i=1:3
            if ~isempty(P1R{i,1}) && ~isempty(P2R{i,1})
                PplusR{i,1} = P1R{i,1} + P2R{i,1};
            elseif isempty(P1R{i,1})
                PplusR{i,1} = P2R{i,1};
            elseif isempty(P2R{i,1})
                PplusR{i,1} = P1R{i,1};
            else
                PplusR{i,1} = [];
            end
        end
        Pplus.(f{:}) = PplusR;
    end
    
    fset = {'Ryy','Ry2','R2y'};
    for f=fset
        P1R = P1.(f{:});    P2R = P2.(f{:});
        PplusR = cell(1,3);
        for j=1:3
            if ~isempty(P1R{1,j}) && ~isempty(P2R{1,j})
                PplusR{1,j} = P1R{1,j} + P2R{1,j};
            elseif isempty(P1R{1,j})
                PplusR{1,j} = P2R{1,j};
            elseif isempty(P2R{1,j})
                PplusR{1,j} = P1R{1,j};
            else
                PplusR{1,j} = [];
            end
        end
        Pplus.(f{:}) = PplusR;
    end
    
    fset = {'R22'};
    for f=fset
        P1R = P1.(f{:});    P2R = P2.(f{:});
        PplusR = cell(1,3);
        for i=1:3
            for j=1:3
                if ~isempty(P1R{i,j}) && ~isempty(P2R{i,j})
                    PplusR{i,j} = P1R{i,j} + P2R{i,j};
                elseif isempty(P1R{i,j})
                    PplusR{i,j} = P2R{i,j};
                elseif isempty(P2R{i,j})
                    PplusR{i,j} = P1R{i,j};
                else
                    PplusR{i,j} = [];
                end
            end
        end
        Pplus.(f{:}) = PplusR;
    end
        

end
end
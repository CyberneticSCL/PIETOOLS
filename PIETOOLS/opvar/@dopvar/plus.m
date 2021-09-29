function [Pplus] = plus(P1,P2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pplus] = plus(P1,P2) performs addition of two operators P: R^p x L2^q to R^m x L2^n
% 
% INPUT
% P1, P2: dopvar class objects
% 
% OUTPUT
% Pplus: returns P1+P2
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - plus
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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




if (isa(P1,'polynomial')||isa(P1,'double')||isa(P1,'dpvar'))&&all(size(P1)==[1,1])&&all(P2.dim(:,1)==P2.dim(:,2))
    Pplus = P2;
    Pplus.P=Pplus.P+P1*eye(P2.dim(1,2)); Pplus.R.R0 = Pplus.R.R0+P1*eye(P2.dim(2,2));
elseif (isa(P2,'polynomial')||isa(P2,'double')||isa(P2,'double'))&&all(size(P2)==[1,1])&&all(P1.dim(:,1)==P1.dim(:,2))
    Pplus = P1;
    Pplus.P=Pplus.P+P2*eye(P1.dim(1,2)); Pplus.R.R0 = Pplus.R.R0+P2*eye(P1.dim(2,2));
elseif ~isa(P1,'dopvar')&& ~isa(P1,'opvar')
    P2dim = P2.dim;
    if all(P2dim(2,:)==[0,0])&&all(size(P1)==P2dim(1,:))
        Pplus = P2;
        Pplus.P = P1+Pplus.P;
    elseif all(P2dim(1,:)==[0,0])&&all(size(P1)==P2dim(2,:))
        Pplus = P2;
        Pplus.R.R0 = P1+Pplus.R.R0;
    end
elseif ~isa(P2,'dopvar')&&~isa(P2,'opvar')
    P1dim = P1.dim;
    if all(P1dim(2,:)==[0,0])&&all(size(P2)==P1dim(1,:))
        Pplus = P1;
        Pplus.P = P2+Pplus.P;
    elseif all(P1dim(1,:)==[0,0])&&all(size(P2)==P1dim(2,:))
        Pplus = P1;
        Pplus.R.R0 = P2+Pplus.R.R0;
    end
else %both are PI operators
    if any(any(P1.dim~=P2.dim))
        error('Operators dimensions do not match. Check the validity of the addition operation.');
    end
    if any(P1.I~=P2.I)
        error('Operators act on different intervals and cannot be added');
    end
    %Create a holder variable for the resultant operator
    Pplus = P1;
    
    
    fset = {'P', 'Q1', 'Q2'};
    
    for i=fset
        if ~isempty(P1.(i{:})) && ~isempty(P2.(i{:}))
            Pplus.(i{:}) = P1.(i{:}) + P2.(i{:});
        elseif isempty(P1.(i{:}))
            Pplus.(i{:}) = P2.(i{:});
        elseif isempty(P2.(i{:}))
            Pplus.(i{:}) = P1.(i{:});
        else
            Pplus.(i{:}) = [];
        end
    end
    
    fset = {'R0','R1','R2'};
    for i=fset
        if ~isempty(P1.R.(i{:})) && ~isempty(P2.R.(i{:}))
            Pplus.R.(i{:}) = P1.R.(i{:}) + P2.R.(i{:});
        elseif isempty(P1.R.(i{:}))
            Pplus.R.(i{:}) = P2.R.(i{:});
        elseif isempty(P2.R.(i{:}))
            Pplus.R.(i{:}) = P1.R.(i{:});
        else
            Pplus.R.(i{:}) = [];
        end
    end
end
end
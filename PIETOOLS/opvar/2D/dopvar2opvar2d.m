function P = dopvar2opvar2d(D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = dopvar2opvar2d converts a dopvar2d object to a corresponding
% opvar2d object
% 
% INPUT
% D: dopvar2d object (or opvar2d object)
% 
% OUTPUT
% P: opvar2d object equivalent of D
%       - Object has same spatial domain and dimensions as input D
%       - All components of P are polynomial class objects
%       - Decision variables appearing in the components of D will be
%         converted to ordinary polynomial variables
%  
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

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
% Initial coding DJ, MP, SS - 08/05/2021

% Error checking
if ~isa(D,'dopvar2d') && ~isa(D,'opvar2d')
    error('Input must be a dopvar2d (or opvar2d) type object')
end

% Initialize the dopvar2d object with the same spatial domain and
% dimensions as the opvar2d object
P = opvar2d([],D.dim,D.I,D.var1,D.var2);

% Convert each component of D to a polynomial class object
fset = {'R00', 'R0x', 'R0y', 'R02', 'Rx0', 'Rxy', 'Ry0', 'Ryx', 'R20'};
for f=fset
    if ~isempty(D.(f{:}))
        P.(f{:}) = polynomial(D.(f{:}));
    else
        P.(f{:}) = polynomial(zeros(size(D.(f{:}))));
    end
end
fset = {'Rxx','Rx2','R2x'};
for f=fset
    tmp_inR = D.(f{:});
    tmp_outR = cell(3,1);
    for i=1:3
        if ~isempty(tmp_inR{i,1})
            tmp_outR{i,1} = polynomial(tmp_inR{i,1});
        else
            tmp_outR{i,1} = polynomial(zeros(size(tmp_inR{i,1})));
        end
    end
    P.(f{:}) = tmp_outR;
end
fset = {'Ryy','Ry2','R2y'};
for f=fset
    tmp_inR = D.(f{:});
    tmp_outR = cell(1,3);
    for j=1:3
        if ~isempty(tmp_inR{1,j})
            tmp_outR{1,j} = polynomial(tmp_inR{1,j});
        else
            tmp_outR{1,j} = polynomial(zeros(size(tmp_inR{1,j})));
        end
    end
    P.(f{:}) = tmp_outR;
end
fset = {'R22'};
for f=fset
    tmp_inR = D.(f{:});
    tmp_outR = cell(3,3);
    for i=1:3
        for j=1:3
            if ~isempty(tmp_inR{i,j})
                tmp_outR{i,j} = polynomial(tmp_inR{i,j});
            else
                tmp_outR{i,j} = polynomial(zeros(size(tmp_inR{i,j})));
            end
        end
    end
    P.(f{:}) = tmp_outR;
end


end
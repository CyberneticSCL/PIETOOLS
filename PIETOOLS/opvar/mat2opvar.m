function Pop = mat2opvar(Pmat,Pdim,vars,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = MAT2OPVAR(PMAT,PDIM,VARS,DOM) takes a matrix or polynomial, and 
% converts it to an 'opvar' object with specified dimensions, variables,
% and domain.
%
% INPUT
% - Pmat:       mxn array of type 'double', 'polynomial', or 'dpvar'
%               representing a multiplier operator to be converted to
%               'opvar' or 'opvar2d' object. If 'dpvar', the output will be 
%               of type 'dopvar' or 'dopvar2d';
% - Pdim:       2x2 (1D) or 4x2 (2D) array of type 'double', specifying the
%               row and column dimensions of the desired 'opvar' or
%               'opvar2d' object. In particular, if
%                   dim = [m0,n0; m1,n1; m2,n2; m3,n3],
%               then the multiplier operator is expected to map
%                   \R^n0 x L2^n1[x] x L2^n2[y] x L2^n3[x,y]
%                   --> \R^m0 x L2^m1[x] x L2^m2[y] x L2^m3[x,y],
%               where n2=n3=m2=m3=0 in the 1D case;
% - vars:       1x2 (1D) or 2x2 (2D) array of type 'polynomial', specifying
%               the primary spatial variables (vars(:,1)=[x;y]) used by the
%               operator, as well as the dummy variables (vars(:,2)) used
%               for integration. Will defaults to [s1,s1_dum] in 1D, and
%               [s1,s1_dum; s2,s2_dum] in 2D;
% - dom:        1x2 (1D) or 2x2 (2D) array specifying for each spatial
%               variable vars(i,1) (and corresponding dummy variable
%               vars(i,2)) the interval [dom(i,1),dom(i,2)] on which this
%               variable exists. Defaults to [0,1] or [0,1;0,1];
%
% OUTPUT
% - Pop:        'opvar' or 'opvar2d' object representing the same
%               multiplier operator as the input Pmat, but now as a PI
%               operator. If 'Pmat' is of type 'dpvar', output will be of
%               type 'dopvar' or 'dopvar2d.
%
% NOTES
% - Currently, the input matrix 'Pmat' is assumed to act solely as
% multiplier operator, mapping e.g. \R-->\R, \R-->L2, or L2-->\R, but not
% mapping e.g. L2-->R, as that would involve an integral operator. As such,
% the function will throw an error if it finds that one of the components
% must be an integral operator in the output structure. To allow conversion
% from matrix to integral operator as well, there is a toggle at the start
% of the function that can be set to true. Note that in that situation,
% e.g. the component mapping L2-->L2 will still be interpreted as a
% multiplier operator, not an integral operator.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  PIETOOLS Team
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
% DJ, 10/23/2024 - Initial coding

% Set toggle to true to allow parts of the matrix to be interpreted as
% kernels in an integral operator.
allow_mat2int = false;

% % % Process the inputs
% % Extract the input arguments.
if nargin==1
    error('Insufficient inputs, input and output dimension of operator are mandatory.')
elseif nargin<=2
    vars = [];
    dom = [];
elseif nargin<=3
    if isa(vars,'double')
        % Assume the third argument specifies the domain
        dom = vars;
        vars = [];
    else
        dom = [];
    end
elseif nargin==4 && (isa(dom,'polynomial') || iscellstr(dom))
    % Assume the third argument specifies the domain
    tmp = dom;
    dom = vars;
    vars = tmp;
end
% % Check that the matrix is properly specified.
use_dopvar = false;
if isa(Pmat,'dpvar')
    use_dopvar = true;
elseif ~isa(Pmat,'double') && ~isa(Pmat,'polynomial')
    error('Matrix to convert must be an object of type ''double'', ''polynomial'', or ''dpvar''.')
end
% % Check that the dimensions are properly specified
if ~isa(Pdim,'double')
    error('Dimensions of the operator should be specified as nx2 array of type double.')
end
if size(Pdim,2)==1 && size(Pmat,1)==size(Pmat,2)
    % Assume a symmetric operator.
    Pdim = [Pdim,Pdim];
elseif size(Pdim,2)~=2
    error('Dimensions of the operator should be specified as nx2 array, providing row dimensions in the first column, and column dimensions in the second.')
end
if any(sum(Pdim,1)~=size(Pmat))
    error('Total number of rows and columns of the matrix should match the total input and output dimensions of the operator.')
end
% % Check that the variables are properly specified.
% Use the default variables if empty 'vars' is specified.
if ~isempty(vars)
    if iscellstr(vars)
        vars = polynomial(vars);
    elseif ~isa(vars,'polynomial') || ~ispvar(vars)
        error('Spatial variables should be specified as Nx2 array of type polynomial, with each element representing an individual variable.')
    end
    if size(vars,2)~=2
        error('Spatial variables should be specified as Nx2 array, providing primary variables in the first column, and dummy variables for integration in the second.')
    end
end
% % Check that the domain is properly specified.
% Use the default domain if empty 'dom' is specified.
if ~isempty(dom)
    if size(dom,2)~=2
        error('Spatial domain should be specified as Nx2 array, providing lower and upper boundaries of the interval on which each spatial variables exists.')
    end
    if any(dom(:,2)<=dom(:,1))
        error('Lower boundary of specified domain exceeds the upper boundary.')
    end
end
% % Check the dimensionality of the operator: 1D or 2D
if ~isempty(vars)
    N = size(vars,1);
elseif ~isempty(dom) && size(dom,1)>1
    N = size(dom,1);
else
    N = log(size(Pdim,1))/log(2);
    if round(N)~=N
        error('Dimension should be specified as array of size [2^N,2], where N is the dimension of the spatial domain on which the operator acts.')
    end
    if N==0
        % If the operator is just R-->R, represent as 1D PI operator with
        % only R-->R component.
        N = 1;
        Pdim = [Pdim;[0,0]];
    end
end
if N>=3
    error('Conversion to operators in more than 2 spatial dimensions is currently not supported.')
end
% Check that an interval is specified for each variable.
if ~isempty(dom) && size(dom,1)==1
    dom = repmat(dom,[N,1]);
elseif ~isempty(dom) && size(dom,1)~=N
    error('The number of rows in the array specifying the domain should match the number of spatial variables.')
end
% Augment dimension array with zeros to match dimension of opvar or opvar2d
% object.
if size(Pdim,1)>2^N
    error('Dimension should be specified as array of size [2^N,2], where N is the dimension of the spatial domain on which the operator acts.')
elseif size(Pdim,1)<2^N
    Pdim = [Pdim;zeros(2^N-size(Pdim,1),2)];
end



% % % Perform the actual conversion.

% Determine the start and end indices of blocks in the matrix that
% correspond to different parameters in the output operator.
nnr = cumsum([0;Pdim(:,1)]);     nnc = cumsum([0;Pdim(:,2)]);
% Distinguish 1D and 2D case.
if N==1
    % % The operator is 1D --> use 'opvar' class.
    % Initialize the operator
    if use_dopvar
        dopvar Pop;
    else
        opvar Pop;
    end
    Pop.dim = Pdim;
    if ~isempty(vars)
        Pop.var1 = vars(1,1);   Pop.var2 = vars(1,2);
    end
    if ~isempty(dom)
        Pop.I = dom;
    end
    % % Set the parameters.
    % Multiplier operator Pop.P:\R-->\R
    Pop.P = Pmat(nnr(1)+1:nnr(2),nnc(1)+1:nnc(2));
    % Integral operator Pop.Q1:L2-->\R
    if allow_mat2int
        Pop.Q1 = Pmat(nnr(1)+1:nnr(2),nnc(2)+1:nnc(3));
    else
        % Currently not allowed, as we interpret matrix only as multiplier
        % operator
        Q1 = Pmat(nnr(1)+1:nnr(2),nnc(2)+1:nnc(3));
        if isa(Q1,'double')
            Q1 = polynomial(Q1);
        end
        if ~isempty(Q1) && any(any(Q1.C))
            error('Conversion of matrix or polynomial to integral operator is currently prohibited.')
        end
    end
    % Multiplier operator Pop.Q2:R-->L2
    Pop.Q2 = Pmat(nnr(2)+1:nnr(3),nnc(1)+1:nnc(2));
    % Multiplier operator Pop.R.R0:L2-->L2
    Pop.R.R0 = Pmat(nnr(2)+1:nnr(3),nnc(2)+1:nnc(3));

else
    % % The operator is 2D --> use 'opvar2d' class.
    % Initialize the operator
    if use_dopvar
        dopvar2d Pop;
    else
        opvar2d Pop;
    end
    Pop.dim = Pdim;
    if ~isempty(vars)
        Pop.var1 = vars(:,1);   Pop.var2 = vars(:,2);
    end
    if ~isempty(dom)
        Pop.I = dom;
    end
    % % Loop over all parameters, setting the parameter equal to the
    % % corresponding block in Pmat.
    params = {'R00', 'R0x', 'R0y', 'R02';
              'Rx0', 'Rxx', 'Rxy', 'Rx2';
              'Ry0', 'Ryx', 'Ryy', 'Ry2';
              'R20', 'R2x', 'R2y', 'R22'};
    % For each parameter params{ii,jj}, isvar(ii,1) indicates whether this
    % parameter maps to a function depending on x, and isvar(ii,2) whether
    % it maps to a function depending on y.
    isvar = [0,0; 1,0; 0,1; 1,1];
    isempty_param = nnr(2:end)==nnr(1:end-1) | (nnc(2:end)==nnc(1:end-1))';
    param_idcs = 1:numel(params);
    for kk=param_idcs(~(isempty_param(:)'))
        [ii,jj] = ind2sub([4,4],kk);
        Pij = Pmat(nnr(ii)+1:nnr(ii+1),nnc(jj)+1:nnc(jj+1));
        if isa(Pij,'double')
            Pij = polynomial(Pij);
        end
        % Check if the parameter maps a function of some variable to a
        % function that does not depend on this variable, i.e., if the
        % parameter defines an integral operator.
        if any(~isvar(ii,:)&isvar(jj,:))
            if ~allow_mat2int && any(any(Pij.C))
                error('Conversion of matrix or polynomial to integral operator is currently prohibited.')
            end
        end
        % Set the actual parameter.
        if isa(Pop.(params{kk}),'cell')
            Pop.(params{kk}){1} = Pij;
        else
            Pop.(params{kk}) = Pij;
        end
    end
end
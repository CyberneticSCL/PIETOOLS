function [Pcat] = vertcat_legacy(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = vertcat(varargin) takes n-inputs and concatentates them vertically,
% provided they satisfy the following criterias.
% 1) Atleast one input is an opvar2d variable.
% 2) If all the inputs are not opvar2d, then the operator maps from R to
% RxL2 or L2 to L2. 
% 3) Currently, it supports RxL2 to RxL2 concatenation only if ALL the inputs are
% opvar2d.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
    Pcat = varargin{1};
else
    a = varargin{1};
    b = varargin{2};

    if isa(a,'opvar2d') % components must have consistent dimensions
        a.dim = a.dim;
    end
    if isa(b,'opvar2d')
        b.dim = b.dim;
    end

    if ~isa(a,'opvar2d') 
        if ~isa(b,'opvar2d')
            if size(a,2)~=size(b,2)
                error('Cannot concatentate vertically. A and B have different input dimensions');
            end
            Pcat = [a;b];
        else
            bdim = b.dim;                
            if size(a,2)~=sum(bdim(:,2))
                error("Cannot concatentate vertically. A and B have different input dimensions");
            end
            Pcat = b;
            Pcat.dim = [b.dim(:,1)+(b.dim(:,1)~=0).*size(a,1),b.dim(:,2)];
            if ~isempty(a)
                if all(bdim(2:4,2) == 0) % a() is from R to R
                    Pcat.R00 = [a; b.R00]; 
                    Pcat.R0x = [zeros(size(a,1),b.dim(2,2)); b.R0x];
                    Pcat.R0y = [zeros(size(a,1),b.dim(3,2)); b.R0y];
                    Pcat.R02 = [zeros(size(a,1),b.dim(4,2)); b.R02];
                elseif all(bdim([1,3:4],2) == 0) % a() is from L2[s1] to L2[s1], Note: For L2 to R, a must be opvar2d 
                    Pcat.Rx0 = [zeros(size(a,1),b.dim(1,2)); b.Rx0];
                    Pcat.Rxx{1,1} = [a; b.Rxx{1,1}];
                    Pcat.Rxx{2,1} = [zeros(size(a,1),b.dim(2,2)); b.Rxx{2,1}];
                    Pcat.Rxx{3,1} = [zeros(size(a,1),b.dim(2,2)); b.Rxx{3,1}];
                    Pcat.Rxy = [zeros(size(a,1),b.dim(3,2)); b.Rxy];
                    Pcat.Rx2{1,1} = [zeros(size(a,1),b.dim(4,2)); b.Rx2{1,1}];
                    Pcat.Rx2{2,1} = [zeros(size(a,1),b.dim(4,2)); b.Rx2{2,1}];
                    Pcat.Rx2{3,1} = [zeros(size(a,1),b.dim(4,2)); b.Rx2{3,1}];
                elseif all(bdim([1:2,4],2) == 0) % a() is from L2[s2] to L2[s2], Note: For L2 to R, a must be opvar2d 
                    Pcat.Ry0 = [zeros(size(a,1),b.dim(1,2)); b.Ry0];
                    Pcat.Ryx = [zeros(size(a,1),b.dim(2,2)); b.Ryx];
                    Pcat.Ryy{1,1} = [a; b.Ryy{1,1}];
                    Pcat.Ryy{1,2} = [zeros(size(a,1),b.dim(3,2)); b.Ryy{1,2}];
                    Pcat.Ryy{1,3} = [zeros(size(a,1),b.dim(3,2)); b.Ryy{1,3}];
                    Pcat.Ry2{1,1} = [zeros(size(a,1),b.dim(4,2)); b.Ry2{1,1}];
                    Pcat.Ry2{1,2} = [zeros(size(a,1),b.dim(4,2)); b.Ry2{1,2}];
                    Pcat.Ry2{1,3} = [zeros(size(a,1),b.dim(4,2)); b.Ry2{1,3}];
                elseif all(bdim(1:3,2) == 0) % a() is from L2[s1,s2] to L2[s1,s2], Note: For L2 to R, a must be opvar2d 
                    Pcat.R20 = [zeros(size(a,1),b.dim(1,2)); b.R20];
                    Pcat.R2x{1,1} = [zeros(size(a,1),b.dim(2,2)); b.R2x{1,1}];
                    Pcat.R2x{2,1} = [zeros(size(a,1),b.dim(2,2)); b.R2x{2,1}];
                    Pcat.R2x{3,1} = [zeros(size(a,1),b.dim(2,2)); b.R2x{3,1}];
                    Pcat.R2y{1,1} = [zeros(size(a,1),b.dim(3,2)); b.R2y{1,1}];
                    Pcat.R2y{1,2} = [zeros(size(a,1),b.dim(3,2)); b.R2y{1,2}];
                    Pcat.R2y{1,3} = [zeros(size(a,1),b.dim(3,2)); b.R2y{1,3}];
                    Pcat.R22{1,1} = [a; b.R22{1,1}];
                    Pcat.R22{2,1} = [zeros(size(a,1),b.dim(4,2)); b.R22{2,1}];
                    Pcat.R22{3,1} = [zeros(size(a,1),b.dim(4,2)); b.R22{3,1}];
                    Pcat.R22{1,2} = [zeros(size(a,1),b.dim(4,2)); b.R22{1,2}];
                    Pcat.R22{2,2} = [zeros(size(a,1),b.dim(4,2)); b.R22{2,2}];
                    Pcat.R22{3,2} = [zeros(size(a,1),b.dim(4,2)); b.R22{3,2}];
                    Pcat.R22{1,3} = [zeros(size(a,1),b.dim(4,2)); b.R22{1,3}];
                    Pcat.R22{2,3} = [zeros(size(a,1),b.dim(4,2)); b.R22{2,3}];
                    Pcat.R22{3,3} = [zeros(size(a,1),b.dim(4,2)); b.R22{3,3}];
                else %find if such an operation is valid in any useful scenario and implement it
                    error('Cannot concatenate vertically. This feature is not yet supported.');
                end
            end
        end
    elseif ~isa(b,'opvar2d') 
        adim = a.dim;
        if size(b,2)~=sum(adim(:,2))
            error("Cannot concatentate vertically. A and B have different input dimensions");
        end
        Pcat = a;
        Pcat.dim = [a.dim(:,1)+(a.dim(:,1)~=0).*size(b,1),a.dim(:,2)];
        if ~isempty(b)
            if all(adim(2:4,2) == 0) % b() is from R to R
                Pcat.R00 = [a.R00; b]; 
                Pcat.R0x = [a.R0x; zeros(size(b,1),a.dim(2,2))];
                Pcat.R0y = [a.R0y; zeros(size(b,1),a.dim(3,2))];
                Pcat.R02 = [a.R02; zeros(size(b,1),a.dim(4,2))];
            elseif all(adim([1,3:4],2) == 0) % b() is from L2[s1] to L2[s1]
                Pcat.Rx0 = [a.Rx0; zeros(size(b,1),a.dim(1,2))];
                Pcat.Rxx{1,1} = [a.Rxx{1,1}; b];
                Pcat.Rxx{2,1} = [a.Rxx{2,1}; zeros(size(b,1),a.dim(2,2))];
                Pcat.Rxx{3,1} = [a.Rxx{3,1}; zeros(size(b,1),a.dim(2,2))];
                Pcat.Rxy = [a.Rxy; zeros(size(b,1),a.dim(3,2))];
                Pcat.Rx2{1,1} = [a.Rx2{1,1}; zeros(size(b,1),a.dim(4,2))];
                Pcat.Rx2{2,1} = [a.Rx2{2,1}; zeros(size(b,1),a.dim(4,2))];
                Pcat.Rx2{3,1} = [a.Rx2{3,1}; zeros(size(b,1),a.dim(4,2))];
            elseif all(adim([1:2,4],2) == 0) % b() is from L2[s2] to L2[s2]
                Pcat.Ry0 = [a.Ry0; zeros(size(b,1),a.dim(1,2))];
                Pcat.Ryx = [a.Ryx; zeros(size(b,1),a.dim(2,2))];
                Pcat.Ryy{1,1} = [a.Ryy{1,1}; b];
                Pcat.Ryy{1,2} = [a.Ryy{1,2}; zeros(size(b,1),a.dim(3,2))];
                Pcat.Ryy{1,3} = [a.Ryy{1,3}; zeros(size(b,1),a.dim(3,2))];
                Pcat.Ry2{1,1} = [a.Ry2{1,1}; zeros(size(b,1),a.dim(4,2))];
                Pcat.Ry2{1,2} = [a.Ry2{1,2}; zeros(size(b,1),a.dim(4,2))];
                Pcat.Ry2{1,3} = [a.Ry2{1,3}; zeros(size(b,1),a.dim(4,2))];
            elseif all(adim(1:3,2) == 0) % b() is from L2[s1,s2] to L2[s1,s2]
                Pcat.R20 = [a.R20; zeros(size(b,1),a.dim(1,2))];
                Pcat.R2x{1,1} = [a.R2x{1,1}; zeros(size(b,1),a.dim(2,2))];
                Pcat.R2x{2,1} = [a.R2x{2,1}; zeros(size(b,1),a.dim(2,2))];
                Pcat.R2x{3,1} = [a.R2x{3,1}; zeros(size(b,1),a.dim(2,2))];
                Pcat.R2y{1,1} = [a.R2y{1,1}; zeros(size(b,1),a.dim(3,2))];
                Pcat.R2y{1,2} = [a.R2y{1,2}; zeros(size(b,1),a.dim(3,2))];
                Pcat.R2y{1,3} = [a.R2y{1,3}; zeros(size(b,1),a.dim(3,2))];
                Pcat.R22{1,1} = [a.R22{1,1}; b];
                Pcat.R22{2,1} = [a.R22{2,1}; zeros(size(b,1),a.dim(4,2))];
                Pcat.R22{3,1} = [a.R22{3,1}; zeros(size(b,1),a.dim(4,2))];
                Pcat.R22{1,2} = [a.R22{1,2}; zeros(size(b,1),a.dim(4,2))];
                Pcat.R22{2,2} = [a.R22{2,2}; zeros(size(b,1),a.dim(4,2))];
                Pcat.R22{3,2} = [a.R22{3,2}; zeros(size(b,1),a.dim(4,2))];
                Pcat.R22{1,3} = [a.R22{1,3}; zeros(size(b,1),a.dim(4,2))];
                Pcat.R22{2,3} = [a.R22{2,3}; zeros(size(b,1),a.dim(4,2))];
                Pcat.R22{3,3} = [a.R22{3,3}; zeros(size(b,1),a.dim(4,2))];
            else %find if such an operation is valid in any useful scenario and implement it
                error('Cannot concatenate vertically. This feature is not yet supported.');
            end
        end
    else
        if any(b.dim(:,2)~=a.dim(:,2))
            error('Cannot concatentate horizontally. A and B have different column dimensions');
        elseif any(any(b.I~=a.I))
            error('Cannot concatentate horizontally: A and B have different domains');
        end
        % Initialize the concatenated operator
        newdim = [a.dim(:,1)+b.dim(:,1),a.dim(:,2)];
        Pcat = opvar2d([],newdim,a.I,a.var1,a.var2);

        % Only concatenate columns which are nonempty
        fset = {};
        c = zeros(4,1);
        if Pcat.dim(1,2)~=0
            fset = [fset,'R00','Rx0','Ry0','R20'];
            c(1) = 1;
        end
        if Pcat.dim(2,2)~=0
            fset = [fset,'R0x','Ryx'];
            c(2) = 1;
        end
        if Pcat.dim(3,2)~=0
            fset = [fset,'R0y','Rxy'];
            c(3) = 1;
        end
        if Pcat.dim(4,2)~=0
            fset = [fset,'R02'];
            c(4) = 1;
        end

        % Perform the concatenation
        for f=fset
            Pcat.(f{:}) = [a.(f{:}); b.(f{:})];
        end
        for i=1:3
            if c(2)
                Pcat.Rxx{i,1} = [a.Rxx{i,1}; b.Rxx{i,1}];
                Pcat.R2x{i,1} = [a.R2x{i,1}; b.R2x{i,1}];
            end
            if c(3)
                Pcat.Ryy{1,i} = [a.Ryy{1,i}; b.Ryy{1,i}];
                Pcat.R2y{1,i} = [a.R2y{1,i}; b.R2y{1,i}];
            end
            if c(4)
                Pcat.Rx2{i,1} = [a.Rx2{i,1}; b.Rx2{i,1}];
                Pcat.Ry2{1,i} = [a.Ry2{1,i}; b.Ry2{1,i}];
                for j=1:3
                    Pcat.R22{i,j} = [a.R22{i,j}; b.R22{i,j}];
                end
            end
        end
    end
    if nargin>2 % Continue concatenation if inputs are more than 2
        Pcat = vertcat(Pcat, varargin{3:end});
    end
end
end
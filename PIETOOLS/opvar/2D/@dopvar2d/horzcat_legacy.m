function [Pcat] = horzcat_legacy(varargin) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = horzcat(varargin) takes n-inputs and concatentates them horizontally,
% provided they satisfy the following criterias.
% 1) Atleast one input is a dopvar2d variable.
% 2) If all the inputs are not dopvar2d, then the operator maps from RxL2 to
%    L2 or R to R. 
% 3) Currently, it supports RxL2 to RxL2 concatenation only if ALL the inputs are
%    dopvar2d.
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

    if isa(a,'dopvar2d') || isa(a,'opvar2d')
        a.dim = a.dim;
    end
    if isa(b,'dopvar2d') || isa(b,'opvar2d')
        b.dim = b.dim;
    end

    if ~isa(a,'dopvar2d') && ~isa(a,'opvar2d')
        if ~isa(b,'dopvar2d') && ~isa(b,'opvar2d') % both are not dopvar2d variables
            if size(a,1)~=size(b,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end
            Pcat = [a b];
        elseif all(b.dim(2:4,1)==0) % a() is from R to R, Note: For L2 to R, a needs to be an opvar2d
            if size(a,1)~=b.dim(1,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end 
            Pcat = b;
            Pcat.R00 = [a b.R00];
            Pcat.Rx0 = [zeros(b.dim(2,1),size(a,2)) b.Rx0];
            Pcat.Ry0 = [zeros(b.dim(3,1),size(a,2)) b.Ry0];
            Pcat.R20 = [zeros(b.dim(4,1),size(a,2)) b.R20];
        elseif all(b.dim([1,3:4],1)==0) % a() is from R to L2[s1], Note: For L2 to L2, a needs to be an opvar2d
            if size(a,1)~=b.dim(2,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end 
            Pcat = b;
            Pcat.R00 = [zeros(b.dim(1,1),size(a,2)) b.R00];
            Pcat.Rx0 = [a b.Rx0];
            Pcat.Ry0 = [zeros(b.dim(3,1),size(a,2)) b.Ry0];
            Pcat.R20 = [zeros(b.dim(4,1),size(a,2)) b.R20];
        elseif all(b.dim([1:2,4],1)==0) % a() is from R to L2[s2], Note: For L2 to L2, a needs to be an opvar2d
            if size(a,1)~=b.dim(3,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end 
            Pcat = b;
            Pcat.R00 = [zeros(b.dim(1,1),size(a,2)) b.R00];
            Pcat.Rx0 = [zeros(b.dim(2,1),size(a,2)) b.Ry0];
            Pcat.Ry0 = [a b.Ry0];
            Pcat.R20 = [zeros(b.dim(4,1),size(a,2)) b.R20];
        elseif all(b.dim(1:3,1)==0) % a() is from R to L2[s1,s2], Note: For L2 to L2, a needs to be an opvar2d
            if size(a,1)~=b.dim(4,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end 
            Pcat = b;
            Pcat.R00 = [zeros(b.dim(1,1),size(a,2)) b.R00];
            Pcat.Rx0 = [zeros(b.dim(2,1),size(a,2)) b.Rx0];
            Pcat.Ry0 = [zeros(b.dim(3,1),size(a,2)) b.Ry0];
            Pcat.R20 = [a b.R20];
        else %find if such an operation is valid in any useful scenario and implement it
            error('Cannot concatenate horizontally. This feature is not yet supported.');
        end
    elseif ~isa(b,'dopvar2d') && ~isa(b,'opvar2d')
        if all(a.dim(2:4,1)==0) % a() is from R to R, Note: L2 to L2 needs to be an opvar
            if size(b,1)~=a.dim(1,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end
            Pcat = a;
            Pcat.R00 = [a.R00 b];
            Pcat.Rx0 = [a.Rx0 zeros(a.dim(2,1),size(b,2))];
            Pcat.Ry0 = [a.Ry0 zeros(a.dim(3,1),size(b,2))];
            Pcat.R20 = [a.R20 zeros(a.dim(4,1),size(b,2))];
        elseif all(a.dim([1,3:4],1)==0) % a() is from R to L2[s1], Note: L2 to L2 needs to be an opvar
            if size(b,1)~=a.dim(2,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end
            Pcat = a;
            Pcat.R00 = [a.R00 zeros(a.dim(1,1),size(b,2))];
            Pcat.Rx0 = [a.Rx0 b];
            Pcat.Ry0 = [a.Ry0 zeros(a.dim(3,1),size(b,2))];
            Pcat.R20 = [a.R20 zeros(a.dim(4,1),size(b,2))];
        elseif all(a.dim([1:2,4],1)==0) % a() is from R to L2[s2], Note: L2 to L2 needs to be an opvar
            if size(b,1)~=a.dim(3,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end
            Pcat = a;
            Pcat.R00 = [a.R00 zeros(a.dim(1,1),size(b,2))];
            Pcat.Rx0 = [a.Rx0 zeros(a.dim(2,1),size(b,2))];
            Pcat.Ry0 = [a.Ry0 b];
            Pcat.R20 = [a.R20 zeros(a.dim(4,1),size(b,2))];
        elseif all(a.dim(1:3,1)==0) % a() is from R to L2[s1,s2], Note: L2 to L2 needs to be an opvar
            if size(b,1)~=a.dim(4,1) 
                error('Cannot concatentate horizontally. A and B have different output dimensions');
            end
            Pcat = a;
            Pcat.R00 = [a.R00 zeros(a.dim(1,1),size(b,2))];
            Pcat.Rx0 = [a.Rx0 zeros(a.dim(2,1),size(b,2))];
            Pcat.Ry0 = [a.Ry0 zeros(a.dim(3,1),size(b,2))];
            Pcat.R20 = [a.R20 b];
        else %find if such a operation is valid is any useful scenario and implement it
            error('Cannot concatenate horizontally.This feature is not yet supported.');
        end
    else % Both arguments are opvar2d or dopvar2d
        if any(b.dim(:,1)~=a.dim(:,1))
            error('Cannot concatentate horizontally. A and B have different row dimensions');
        elseif any(any(b.I~=a.I))
            error('Cannot concatentate horizontally: A and B have different domains');
        end
        % Initialize the concatenated operator
        newdim = [a.dim(:,1),a.dim(:,2)+b.dim(:,2)];
        Pcat = dopvar2d([],newdim,a.I,a.var1,a.var2);

        % Only concatenate rows which are nonempty
        fset = {};
        r = zeros(4,1);
        if Pcat.dim(1,1)~=0
            fset = [fset,'R00','R0x','R0y','R02'];
            r(1) = 1;
        end
        if Pcat.dim(2,1)~=0
            fset = [fset,'Rx0','Rxy'];
            r(2) = 1;
        end
        if Pcat.dim(3,1)~=0
            fset = [fset,'Ry0','Ryx'];
            r(3) = 1;
        end
        if Pcat.dim(4,1)~=0
            fset = [fset,'R20'];
            r(4) = 1;
        end

        % Perform the concatenation
        for f=fset
            Pcat.(f{:}) = [a.(f{:}) b.(f{:})];
        end
        for i=1:3
            if r(2)
                Pcat.Rxx{i,1} = [a.Rxx{i,1} b.Rxx{i,1}];
                Pcat.Rx2{i,1} = [a.Rx2{i,1} b.Rx2{i,1}];
            end
            if r(3)
                Pcat.Ryy{1,i} = [a.Ryy{1,i} b.Ryy{1,i}];
                Pcat.Ry2{1,i} = [a.Ry2{1,i} b.Ry2{1,i}];
            end
            if r(4)
                Pcat.R2x{i,1} = [a.R2x{i,1} b.R2x{i,1}];
                Pcat.R2y{1,i} = [a.R2y{1,i} b.R2y{1,i}];
                for j=1:3
                    Pcat.R22{i,j} = [a.R22{i,j} b.R22{i,j}];
                end
            end
        end
    end
    if nargin>2 % check if there are more than 2 objects that need to stacked
        Pcat = horzcat(Pcat, varargin{3:end});
    end
end
end
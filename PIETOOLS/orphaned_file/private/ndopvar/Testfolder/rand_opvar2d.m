function Pop = rand_opvar2d(dim,deg,dom,var1,var2)
% RAND_OPVAR2D Generates a random 2D-PI operator with monomials of degree 
% at most deg in each of the variables

% Initialize an empty operator
Pop = opvar2d(); 
% Check that the dimensions are properly specified
if numel(dim)==2
    dim = [0,0;0,0;0,0;dim(1),dim(2)];
elseif ~all(size(dim)==[4,2])
    error("Dimensions should be specified as 4x2 array.")
end
Pop.dim = dim;
m0 = dim(1,1);      n0 = dim(1,2);
mx = dim(2,1);      nx = dim(2,2);
my = dim(3,1);      ny = dim(3,2);
m2 = dim(4,1);      n2 = dim(4,2);

% Set the spatial domain
if nargin>=3
    Pop.I = dom;
end
% Set the variables
if nargin>=4
    Pop.var1 = var1;
else
    var1 = Pop.var1;
end
if nargin>=5
    Pop.var2 = var2;
else
    var2 = Pop.var2;
end
s1 = var1(1);       t1 = var2(1);
s2 = var1(2);       t2 = var2(2);
            
% % Set the randomly generated parameters
% Start with parameters consisting of a single term
Pop.R00 = rand([m0,n0]);
Pop.R0x = rand_poly([m0,nx],s1,deg);
Pop.R0y = rand_poly([m0,ny],s2,deg);
Pop.R02 = rand_poly([m0,n2],[s1;s2],deg);
Pop.Rx0 = rand_poly([mx,n0],s1,deg);
Pop.Ry0 = rand_poly([my,n0],s2,deg);
Pop.R20 = rand_poly([m2,n0],[s1;s2],deg);
Pop.Rxy = rand_poly([mx,ny],[s1;s2],deg);
Pop.Ryx = rand_poly([my,nx],[s1;s2],deg);
% Next, set 3-PI operators along x-direction
Pop.Rxx{1} = rand_poly([mx,nx],s1,deg);
Pop.Rx2{1} = rand_poly([mx,n2],[s1;s2],deg);
Pop.R2x{1} = rand_poly([m2,nx],[s1;s2],deg);
for i=2:3
    Pop.Rxx{i} = rand_poly([mx,nx],[s1;t1],deg);
    Pop.Rx2{i} = rand_poly([mx,n2],[s1;s2;t1],deg);
    Pop.R2x{i} = rand_poly([m2,nx],[s1;s2;t1],deg);
end
% Next, set 3-PI operators along y-direction
Pop.Ryy{1} = rand_poly([my,ny],s2,deg);
Pop.Ry2{1} = rand_poly([my,n2],[s1;s2],deg);
Pop.R2y{1} = rand_poly([m2,ny],[s1;s2],deg);
for i=2:3
    Pop.Ryy{i} = rand_poly([my,ny],[s2;t2],deg);
    Pop.Ry2{i} = rand_poly([my,n2],[s1;s2;t2],deg);
    Pop.R2y{i} = rand_poly([m2,ny],[s1;s2;t2],deg);
end
% Finally, set 9-PI operator
Pop.R22{1,1} = rand_poly([m2,n2],[s1;s2],deg);
for ii=2:3
    Pop.R22{ii,1} = rand_poly([m2,n2],[s1;s2;t1],deg);
    Rop.R22{ii,1} = rand_poly([n2,n2],[s1;s2;t1],deg);
    Pop.R22{1,ii} = rand_poly([m2,n2],[s1;s2;t2],deg);
    Rop.R22{1,ii} = rand_poly([n2,n2],[s1;s2;t2],deg);
    for jj=2:3
        Pop.R22{ii,jj} = rand_poly([m2,n2],[s1 s2 t1 t2],deg);
        Rop.R22{ii,jj} = rand_poly([n2,n2],[s1;s2;t1;t2],deg);
    end
end


end


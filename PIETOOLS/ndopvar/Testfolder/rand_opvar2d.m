function Pop = rand_opvar2d(dim,deg,dom,var1,var2)
% RAND_OPVAR2D Generates a random 9-PI operator with monomials of degree at
% most deg in each of the variables

Pop = opvar2d(); 
if numel(dim)~=2
    error("Only 9-PI operators are supported.")
end
m = dim(1);     n = dim(2);
Pop.dim = [0,0;0,0;0,0;m,n]; 

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
s1 = var1(1);       s1_dum = var2(1);
s2 = var1(2);       s2_dum = var2(2);
            
% Set the randomly generated parameters
Pop.R22{1,1} = rand_poly([m,n],[s1;s2],deg);
for ii=2:3
    Pop.R22{ii,1} = rand_poly([m,n],[s1;s2;s1_dum],deg);
    Rop.R22{ii,1} = rand_poly([n,n],[s1;s2;s1_dum],deg);
    Pop.R22{1,ii} = rand_poly([m,n],[s1;s2;s2_dum],deg);
    Rop.R22{1,ii} = rand_poly([n,n],[s1;s2;s2_dum],deg);
    for jj=2:3
        Pop.R22{ii,jj} = rand_poly([m,n],[s1 s2 s1_dum s2_dum],deg);
        Rop.R22{ii,jj} = rand_poly([n,n],[s1;s2;s1_dum;s2_dum],deg);
    end
end


end


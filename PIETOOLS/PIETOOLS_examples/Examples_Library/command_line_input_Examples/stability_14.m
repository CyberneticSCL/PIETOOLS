clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
%   PDE: x1_{t} = (1/aa) * x2_{s}                
%        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3   
%        x3_{t} = (1/II)*x2 + (1/II)*x4_{s}     
%        x4_{t} = II * x3_{s}                   
%   BCs: x1(0) = 0,         x2(1) = 0           
%        x3(0) = 0,         x4(1) = 0           
x1=state('pde');x2=state('pde');x3=state('pde');x4=state('pde'); 
pde = sys(); 
k = 1;    aa = 1;    II = 1;    g = 1;    
%% Define equations
%   PDE: x1_{t} = (1/aa) * x2_{s}                
%        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3   
%        x3_{t} = (1/II)*x2 + (1/II)*x4_{s}     
%        x4_{t} = II * x3_{s}
eq_PDE = [diff(x1,t)==(1/aa)*diff(x2,s); 
          diff(x2,t)==k*aa*g*diff(x1,s)-k*aa*g*x3;
          diff(x3,t)==(1/II)*x2+(1/II)*diff(x4,s);
          diff(x4,t)==II*diff(x3,s)];
eq_BC = [subs(x1,s,0)==0;subs(x3,s,0)==0;
         subs(x2,s,1)==0;subs(x4,s,1)==0]; %   BCs:BCs: x1(0) = 0, x2(1) = 0, x3(0) = 0, x4(1) = 0 
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);

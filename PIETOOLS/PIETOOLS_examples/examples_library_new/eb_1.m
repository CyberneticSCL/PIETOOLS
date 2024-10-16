% % PDE        u_{tt} = -c*u_{ssss}
% % with BCs   u(s=0) = 0
% %            u_{s}(s=0) = 0 
% %            u_{ss}(s=1) = 0 
% %            u_{sss}(s=1) = 0 
function obj = eb_1()
pvar t s1;
pdevar x;
eqns = [diff(x,t,2)==-c*diff(x,s1,4);
             subs(x,s1,0)==0; subs(diff(x,s1),s1,0)==0;
             subs(diff(x,s1,2),s1,1)==0; subs(diff(x,s1,3),s1,1)==0];
obj = sys('pde',eqns);
end
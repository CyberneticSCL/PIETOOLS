n=7;
T = sdpvar(3,n); T2 = sdpvar(n); gam = sdpvar(1);
B = [[-gam,    0,    0;
     0,    -gam,    1;
     0 ,      1,    0], T; T', T2];
optimize(B<=0,[])

%%
n=9;
T = sdpvar(4,n); T2 = sdpvar(n);gam = sdpvar(1);


C = [[-gam,    0,    0,    0; 
     0,    -gam,    1,    1; 
     0 ,      1,    0,    0; 
     0 ,      1,    0,    0], T; T', T2];



optimize(C<=0,[])



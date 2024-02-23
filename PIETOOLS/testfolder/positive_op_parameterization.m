pvar s theta;
T = opvar({0,0;0,{0,0,0}});
Q = T'*T;
n=5;
for i=0:n^2
    eval("pvar a"+num2str(i));
end

opvar Z P;
Z.R.R0 = [1;s;0;0;0]; Z.R.R1 = [0;0;1; s; 0]; Z.R.R2 = [0;0;1; 0; theta];
P.R.R0 = [a0   a1   a2   a3   a4; 
                a5   a6   a7   a8   a9;
                a10 a11 a12 a13 a14;
                a15 a16 a17 a18 a19;
                a20 a21 a22 a23 a24;];
Qd = Z'*P*Z;
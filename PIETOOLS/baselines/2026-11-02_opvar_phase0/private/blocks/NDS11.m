% examples_NDSDDF_library_PIETOOLS.m lines 190-203, one comment level removed (lb_blocks)
 stability=1
 stability_dual=1
 
 ki=10;kp=10;
 a=.4;b=50;h=.2;d=.8;sig=.3;
 al1=d+kp;gam1=b*ki*d+a*ki;
 al2=(d-kp)*exp(sig*h);gam2=(b*ki*d-a*ki)*exp(sig*h);
 beta1=(b*kp+a)*d+b*d^2+a*kp+ki;
 beta2=((b*kp+a)*d-b*d^2-a*kp-ki)*exp(sig*h);
 
 NDS.A0=(1/al1)*[0 al1;-sig^2*al1+sig*beta1-gam1 -beta1+2*sig*al1];
 NDS.Ai{1}=(1/al1)*[0 0;-sig^2*al2+sig*beta2-gam2 -beta2+2*sig*al2];
 NDS.Ei{1}=[0 0;0 -al2/al1];
 NDS.tau(1)=h;

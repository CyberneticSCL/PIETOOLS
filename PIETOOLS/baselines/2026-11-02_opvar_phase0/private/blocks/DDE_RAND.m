% examples_DDE_library_PIETOOLS.m lines 444-498, one comment level removed (lb_blocks)
 seed='default';%234235;
 maxn=3; pvar_deg=1;
 stability=0;
 stability_dual=0;
 Hinf_gain=0;
 Hinf_gain_dual=0;
 Hinf_control=0;
 Hinf_estimator=1;
 
 %%%%%%%% Begin generating case
 nK=randi(maxn+1)-1; nx =randi(maxn+1)-1; nw=randi(maxn+1)-1; nu=randi(maxn+1)-1; nz=randi(maxn+1)-1; ny=randi(maxn+1)-1;
 %nK=2; nx =2; nw=1; nu=2; nz=1; ny=2;
 rng(seed)
 DDE.tau=rand([1 nK]);
 DDE.A0=rand([nx nx]);
 DDE.B1=rand([nx nw]);
 DDE.B2=rand([nx nu]);
 DDE.C1=rand([nz nx]);
 DDE.C2=rand([ny nx]);
 DDE.D11=rand([nz nw]);
 DDE.D12=rand([nz nu]);
 DDE.D21=rand([ny nw]);
 DDE.D22=rand([ny nu]);
 for i=1:nK
     DDE.Ai{i}=rand([nx nx]);
     DDE.B1i{i}=rand([nx nw]);
     DDE.B2i{i}=rand([nx nu]);
     DDE.C1i{i}=rand([nz nx]);
     DDE.C2i{i}=rand([ny nx]);
     DDE.D11i{i}=rand([nz nw]);
     DDE.D12i{i}=rand([nz nu]);
     DDE.D21i{i}=rand([ny nw]);
     DDE.D22i{i}=rand([ny nu]);
     DDE.Adi{i}=rand([nx nx]);
     DDE.B1di{i}=rand([nx nw]);
     DDE.B2di{i}=rand([nx nu]);
     DDE.C1di{i}=rand([nz nx]);
     DDE.C2di{i}=rand([ny nx]);
     DDE.D11di{i}=rand([nz nw]);
     DDE.D12di{i}=rand([nz nu]);
     DDE.D21di{i}=rand([ny nw]);
     DDE.D22di{i}=rand([ny nu]);
     for j=1:pvar_deg
         DDE.Adi{i}=DDE.Adi{i}+rand([nx nx])*s^j;
         DDE.B1di{i}=DDE.B1di{i}+rand([nx nw])*s^j;
         DDE.B2di{i}=DDE.B2di{i}+rand([nx nu])*s^j;
         DDE.C1di{i}=DDE.C1di{i}+rand([nz nx])*s^j;
         DDE.C2di{i}=DDE.C2di{i}+rand([ny nx])*s^j;
         DDE.D11di{i}=DDE.D11di{i}+rand([nz nw])*s^j;
         DDE.D12di{i}=DDE.D12di{i}+rand([nz nu])*s^j;
         DDE.D21di{i}=DDE.D21di{i}+rand([ny nw])*s^j;
         DDE.D22di{i}=DDE.D22di{i}+rand([ny nu])*s^j;
     end
     
 end

% examples_DDE_library_PIETOOLS.m lines 220-241, one comment level removed (lb_blocks)
 % % Scalability Example from [4]
 % delaymat=[1 2 3 5 10]
 % dimmat=[1 2 3 5 10]
 %  for iii=5:5 %dimmat=[1 2 3 5 10] n
 %      for jjj=1:5 %delaymat=[1 2 3 5 10] K
 %         clear prog tau Ai A0 B1 B2 C1 C1i D11 D12 D21 D22 Dop
 %  nndim=dimmat(iii);nndelay=delaymat(jjj);
 
 Hinf_control=1
  nndim=3;nndelay=3;
 
 DDE.A0=zeros(nndim);
 DDE.B1=ones(nndim,1);
 DDE.B2=ones(nndim,1);
 DDE.C1=[ones(1,nndim);zeros(1,nndim)];
 DDE.D11=[0;0]; DDE.D12=[0;1];
 for j=1:nndelay
   DDE.Ai{j}=-eye(nndim)/nndelay;
 %  DDE.C1i{j}=zeros(2,nndim);
   DDE.tau(j) = j/nndelay;
 end
% Answer [3,3] - .9630

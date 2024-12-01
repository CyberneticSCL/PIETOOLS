function out = blkdiag(A,B)
dimAin = A.dim(:,2);
dimAout = A.dim(:,1);
dimBin = B.dim(:,2);
dimBout = B.dim(:,1);
zeroAB = buildopvar('dim',[dimAout,dimBin],'dom',B.dom);
zeroBA = buildopvar('dim',[dimBout,dimAin],'dom',A.dom);
out = [A,zeroAB; zeroBA,B];
end
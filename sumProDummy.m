
load('out1_1')
-2*mean(Ls)
mat1 = [x0; mean(matpara,1); quantile(matpara,[.025,.975],1)]';

load('out2_0')
-2*mean(Ls)
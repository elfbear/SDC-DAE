function [x1vec,x2vec,yvec]= exactSol(w,tvec)
x1vec = exp(2*w*tvec');
x2vec = exp(w*tvec');
yvec = w*x2vec-x1vec;
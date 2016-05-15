function [x1vec,x2vec,yvec]= exactSol(w,tvec)
x1vec = exp(w*tvec');
x2vec = exp(w*tvec');
yvec = w*exp(w*tvec')/2;
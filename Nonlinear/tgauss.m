function [A,B,tc]=tgauss(n)
%A=zeros(n+1,n+1);
eval(strcat('tgauss',num2str(n)));
return
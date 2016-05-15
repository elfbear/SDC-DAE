function [vect,A,B] = GaussNodes(t0, t1, p)
addpath 'nodes'
% tc is the gauss nodes between [0,1]
[A B tc]=tgauss(p);
% vect returns gauss nodes between [0,1] and the end pts
if(tc(end)~=1)
    vect = [(t1-t0)*tc+t0,t1];
else
    vect = (t1-t0)*tc+t0;
end

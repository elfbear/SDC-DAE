% This code is for testing the hybrid parallel + SDC algorithm for
% semi-explicit DAE
% By Wenhua Guan @ pnnl, July 22nd, 2015

% example :  x1'(t) = w * x2(t)
%            x2'(t) = w/2* x1(t)+ y(t)
%            y(t)   = w*x1(t)/2

clc, clear all, close all
w = .1;
t0 = 0; tfinal = 5; dt = 1;
N = (tfinal-t0)/dt;

% Initials
x1_0=1; x2_0 =1; y0 = w/2;
p = 5;

% coarse grids 
tvec = t0:dt:tfinal;   

%step1: Use trapezoid rule to compute the low-order provisional solution
[x1vec0, x2vec0, yvec0] = LowOrder(w,tvec,x1_0,x2_0,y0);
[x1vec, x2vec, yvec]= exactSol(w,tvec);

error1 = norm(x1vec-x1vec0);
error2 = norm(x2vec-x2vec0);
error3 = norm(yvec -yvec0);
fprintf('Predictor: %d, %d, %d\n', error1, error2, error3);

%step2: Precorrector 
[x1mat0, x2mat0, ymat0, delta1_0, delta2_0] = PreCorrector(p,N,w,tvec,...
                                                x1vec0,x2vec0,yvec0);

%step 3 : Corrector
M= 7;
[x1mat, x2mat, ymat] = Corrector(w,tvec, p, N, M, x1mat0, x2mat0, ymat0);

%Check and plot error
errx1 = abs(x1vec(2:end)-x1mat(:,end));
errx2  = abs(x2vec(2:end)-x2mat(:,end));
erry =  abs(yvec(2:end) - ymat(:,end));
ploterror(tvec, errx1,errx2,erry);






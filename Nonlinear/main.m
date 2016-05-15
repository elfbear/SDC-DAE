% This code is for testing the hybrid parallel + SDC algorithm for
% semi-explicit nonliear DAE
% By Wenhua Guan @ pnnl, July 31, 2015

% nonliear example :  x1'(t) = 2 * w * x2(t)^2
%                     x2'(t) = x1(t)+ y(t)
%                     y(t)   = w*x2(t)- x1(t)

clc, clear all, close all
format long
w = .1;
t0 = 0; tfinal = 5; dt = 1;
N = (tfinal-t0)/dt;

% Initials
x1_0=1; x2_0 =1; y0 = w-1;
tvec = t0:dt:tfinal;   % coarse grids 
% p = 5;
errx1 = zeros(8,1);errx2 = errx1;erry = errx1;
for p = 3:10
%step1: Use trapezoid rule to compute the low-order provisional solution
[x1vec0, x2vec0, yvec0] = LowOrder(w,tvec,x1_0,x2_0,y0);
[x1vec, x2vec, yvec]= exactSol(w,tvec);

error1 = (x1vec(2:end)-x1vec0(2:end));
error2 = (x2vec(2:end)-x2vec0(2:end));
error3 = (yvec(2:end) -yvec0(2:end));
fprintf('Predictor: %d, %d,%d\n', norm(error1,Inf), norm(error2, Inf), ...
       norm(error3, Inf));

% %step2: Precorrector 
[x1mat0, x2mat0, ymat0, delta1_0, delta2_0] = PreCorrector(p,N,w,tvec,...
                                                 x1vec0,x2vec0,yvec0);

%step 3 : Corrector
% errx1 = zeros(13,1);errx2 = errx1;erry = errx1;
% for M= 3:15
M=5;
[x1mat, x2mat, ymat] = Corrector(w,tvec, p, N, M, x1mat0, x2mat0, ymat0);

%Check and plot error
errx1(p-2) = abs(x1vec(end)-x1mat(end,end));
errx2(p-2)  = abs(x2vec(end)-x2mat(end,end));
erry(p-2) =  abs(yvec(end) - ymat(end,end));
end
% ploterror(tvec, errx1,errx2,erry);
%  ploterrorM(3:15, -log10(errx1),-log10(errx2),-log10(erry))
ploterrorp(3:10, -log10(errx1),-log10(errx2),-log10(erry))





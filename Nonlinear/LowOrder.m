function [x1vec0, x2vec0, yvec0] = LowOrder(w, tvec, x1_0,x2_0, y0)
N = length(tvec)-1;
%J is the Jacobi matrix
J = eye(3);
J(3,1) = 1;
J(3,2) =-w;

x1vec0 = zeros(N+1,1);
x2vec0 = zeros(N+1,1);
yvec0 = zeros(N+1,1);
    
x1vec0(1) = x1_0;
x2vec0(1) = x2_0;
yvec0(1) = y0;

for i=1:N
    dt =tvec(i+1)-tvec(i);
    J(1,2)= -2*w*dt; % Just part of J(1,2)
    J(2,1) = -dt/2;
    J(2,3) = -dt/2;
%     fprintf('%d time step', i)
    
    % Newton method for Low Order part
    eps = 1e-16;
    maxiter = 50;
    x10 = x1vec0(i);
    x20 = x2vec0(i);
    y0 =  yvec0(i);
    
    [x1vec0(i+1), x2vec0(i+1), yvec0(i+1)]= LowOrderNewton(w,dt, ...
      x1vec0(i),x2vec0(i), yvec0(i), eps, maxiter,J,x10, x20, y0);
end 
%Newton method for the nonlinear part
%This Newton method has been applied for two parts of the hybrid parareal sdc
%algorithm.
%1. Prediction and Precorrector : When implicit low order method has been used;
%   Here, we applied trapezoid rule in the code. LowOrderNewton()
%2, Corrector : Only if backward Euler method is used in SDC. SDCNewton()

function [x1,x2,y] = LowOrderNewton(w,dt, x1n,x2n, yn, eps, maxiter,...
                                    J,x10, x20, y0)
iter=1;
fold = 10*ones(3,1);
while (iter<=maxiter) 
    %fprintf('%d newton iteration', iter)
    f = F(w,dt,x1n, x2n, yn, x10,x20,y0);
    if(norm(f, Inf)<eps | (norm(f, Inf)>= norm(fold, Inf)) )
        x1 = x10; x2 =x20; y = y0;
        break; 
    end
    J(1,2) = J(1,2)*x20;
    
    % inner loop
    delta = -J\f;
    x1 = x10 + delta(1);
    x2 = x20 + delta(2);
    y = y0 + delta(3);
    
    x10 = x1;
    x20 = x2;
    y0 = y;
    
    iter = iter+1;
    fold = f;
end
%To evaluate the Fvalues at t_n
function f = F(w,dt, x1n, x2n, yn,x1,x2,y)
f = zeros(3,1);
f(1) = x1 - w*dt * x2^2 -(x1n + dt*w *x2n^2) ;
f(2) = -dt/2*x1 + x2 -dt/2*y -(x2n+ dt/2*x1n + dt/2 *yn) ;
f(3) = x1 - w*x2+ y;
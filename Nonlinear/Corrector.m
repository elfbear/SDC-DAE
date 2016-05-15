%Corrector :
function [x1mat, x2mat, ymat] = Corrector(w,tvec,p, N, M,...
                                          x1mat0, x2mat0, ymat0)
for m =1:M
    x1mat = zeros(N,p+2);
    x2mat = zeros(N,p+2);
    ymat = zeros(N,p+2);
    
    %communication
    x1mat(1,1) = x1mat0(1,1);
    x2mat(1,1) = x2mat0(1,1);
    ymat(1,1) = ymat0(1,1);

    x1mat(2:end,1) = x1mat0(1:end-1,end);
    x2mat(2:end,1) = x2mat0(1:end-1,end);
    ymat(2:end,1) = ymat0(1:end-1,end);
%     fprintf('start %d sdc iteration\n',m)
    
    % To compute [t0, tfinal] simultaneously 
    for n = 1 : N 
        % SDC for [t_n, t_{n+1}]
        [x1mat(n,:), x2mat(n,:), ymat(n,:)] = sdc(w,tvec(n),...
            tvec(n+1), p, x1mat(n,1),x2mat(n,1),ymat(n,1),...
            x1mat0(n,:), x2mat0(n,:), ymat0(n,:));  
    end
    x1mat0 = x1mat;
    x2mat0 = x2mat;
    ymat0 = ymat;

%     fprintf('End of %d sdc iteration\n',m)
%     x1vec(2:end)- x1mat(:,end)
%     x2vec(2:end) - x2mat(:,end)
%     yvec(2:end) - ymat(:,end)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SDC procedure for [t0, t1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x1vec,x2vec,yvec] = sdc(w,t0, t1, p,x10,x20,y0,x1vec0, x2vec0, yvec0)
x1vec = zeros(p+2,1);
x2vec = zeros(p+2,1);
yvec  = zeros(p+2,1);
x1vec(1) = x10;
x2vec(1) = x20;
yvec(1) = y0;
% compute the integration matrix
[tnGaussVec,eps1,eps2,B] = epsilion(w,t0,t1,p, x1vec0, x2vec0, yvec0); 

for j =1:p
    dt = tnGaussVec(j+1)- tnGaussVec(j);
    % For non-stiff system: forward Euler. 
%     [x1vec(j+1), x2vec(j+1), yvec(j+1)] = forward_Euler(w,dt, x1vec0(j),...
%          x2vec0(j), yvec0(j), x1vec(j), x2vec(j),yvec(j), eps1(j),eps2(j));
    %For stiff system: backward Euler
    [x1vec(j+1), x2vec(j+1), yvec(j+1)] = backward_Euler(w,dt,x1vec0(j+1),...
         x2vec0(j+1), yvec0(j+1), x1vec(j), x2vec(j),yvec(j), eps1(j),eps2(j));
end
% interpolation on t_{n+1} using the p gauss quadrature pts in [tn, t_{n+1}]
x1vec(p+2) = (t1-t0)* B*x1vec(2:p+1); 
x2vec(p+2) = (t1-t0)* B*x2vec(2:p+1);
yvec(p+2) = gval(w,x1vec(p+2),x2vec(p+2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Euler or backward Euler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For non-stiff system: forward Euler.
function [x1new, x2new, ynew] = forward_Euler(w,dt,x10,x20,y0,...
                                              x1, x2,y,eps1,eps2)
    f1 = fva1(w,x2);
    f2 = fva1(w,x20);
    f3 = fva2(x1, y);
    f4 = fva2(x10,y0);
    x1new = x1 + dt*(f1-f2) + eps1;
    x2new = x2 + dt*(f3-f4) + eps2;
    ynew = gval(w,x1new,x2new);
% For stiff system: backward Euler 
function [x1new, x2new, ynew] = backward_Euler(w,dt,x10,x20,y0,...
                                              x1, x2,y,eps1,eps2)
J = eye(3,3); % The Jacobi matrix of SDC procedure
J(1,2) = -4*w*dt; 
J(2,1) = -dt;
J(2,3) = -dt;
J(3,2) = -w;
J(3,1) = 1;
% parameters for newton method in SDC,use x10,x20,y0 as the initials
tol = 1e-16;
maxiter = 50;
[x1new, x2new, ynew] =  SDCNewton(w,dt, tol,maxiter,x1,x2,x10,x20,y0,...
                                 J,eps1,eps2);
                             
%Newton method for the nonlinear part
function [x1,x2,y] = SDCNewton(w,dt,tol, maxiter, x1j,x2j,x1m,x2m,ym,... 
                                    J,g1,g2)
iter=1;
x10 = x1m;
x20 = x2m;
y0 = ym;
fold = 10*ones(3,1);
while (iter<=maxiter) 
    %fprintf('%d newton iteration', iter)
    f = Fsdc(w,dt,x1j,x2j,x1m, x2m, ym,x10,x20,y0,g1,g2);
    if(norm(f, Inf)<tol | (norm(f, Inf)>= norm(fold, Inf)) )
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

function f = Fsdc(w,dt,x1j,x2j,x1m, x2m, ym,x1,x2,y,g1,g2)
f = zeros(3,1);
f(1) = x1-2*w *dt *x2^2-(x1j - 2*w*dt*x2m^2 + g1);
f(2) = -dt*x1 + x2 -dt*y - (x2j-dt*x1m-dt*ym + g2);
f(3) = x1-w*x2+y;

function [tnGaussVec,eps1, eps2, B] = epsilion(w,t0, t1, p, ...
                                               x1vec0, x2vec0, yvec0)
[tnGaussVec,A,B] = GaussNodes(t0, t1, p);
f1 = fva1(w,x2vec0(2:p+1));
f2 = fva2(x1vec0(2:p+1), yvec0(2:p+1));
eps1 = (t1-t0) * A*f1;
eps2 = (t1-t0) * A*f2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To evaluate the right hand side of x1' = f1(x2) = 2*w*x_2^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f1 = fva1(w,x2)
% f1= zeros(length(x2),1);
f1 = (2* w* x2.^2)';
% To evaluate the right hand side of x2' = f2(x1,y) = x1+y
function f2 = fva2(x1,y)
% f2 = zeros(length(y),1);
f2 = x1' + y';
% To evaluate y = g(x1,x2) = w*x2-x1
function gval= gval(w,x1,x2)
gval = w*x2 -x1;

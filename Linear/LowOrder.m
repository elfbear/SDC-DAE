function [x1vec0, x2vec0, yvec0] = LowOrder(w,tvec, x1_0,x2_0, y0)
N = length(tvec)-1;
A = zeros(3);
A(1,1)=1;
A(3,3)=1;
A(3,1) = -w/2;

x1vec0 = zeros(N+1,1);
x2vec0 = zeros(N+1,1);
yvec0 = zeros(N+1,1);
    
x1vec0(1) = x1_0;
x2vec0(1) = x2_0;
yvec0(1) = y0;

b = zeros(3,1);

for i=1:N
    dt =tvec(i+1)-tvec(i);
    A(1,2)= -w*dt/2;
    A(2,2) = 1-w*dt/4;
    A(2,3) = -dt/2;

    b(1) = x1vec0(i)+w*dt*x2vec0(i)/2 ;
    b(2) = (1+w*dt/4)*x2vec0(i)+dt*yvec0(i)/2;
        
    sol = zeros(3,1); 
    sol = A\b;
    
    x1vec0(i+1) = sol(1);
    x2vec0(i+1) = sol(2);
    yvec0(i+1) = sol(3);
end
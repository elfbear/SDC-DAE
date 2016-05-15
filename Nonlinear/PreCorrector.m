function [x1mat0, x2mat0, ymat0, delta1_0, delta2_0] = PreCorrector(p, N,...
                            w,tvec, x1vec0,x2vec0,yvec0)
x1mat0 = zeros(N, p+2);
x2mat0 = zeros(N,p+2);
ymat0 = zeros(N,p+2);

x1mat0(:,1) = x1vec0(1:N);
x2mat0(:,1) = x2vec0(1:N);
ymat0(:,1) = yvec0(1:N);

err1 = zeros(N,1);
err2 = zeros(N,1);
err3 = zeros(N,1);

% Can Use more processors here
for n = 0: N-1
%     x1vec = zeros(p+2,1); x2vec = zeros(p+2,1); yvec = zeros(p+2,1);
    [tnGaussVec,A,B] = GaussNodes(tvec(n+1), tvec(n+2), p);
    [x1mat0(n+1,:), x2mat0(n+1,:), ymat0(n+1,:)] = LowOrder(w,tnGaussVec,...
                              x1mat0(n+1,1),x2mat0(n+1,1),ymat0(n+1,1));
    [x1vec,x2vec,yvec]= exactSol(w,tnGaussVec);
%     fprintf('%dth time interval\n',n)
%     x1mat0(n+1,:)'-x1vec
%     x2mat0(n+1,:)'-x2vec
%     ymat0(n+1,:)'-yvec
    
    err1(n+1)= norm(x1mat0(n+1,:)'-x1vec, Inf);
    err2(n+1)= norm(x2mat0(n+1,:)'-x2vec, Inf);
    err3(n+1)= norm(ymat0(n+1,:)'-yvec, Inf);
end

maxerr1 = max(err1);
maxerr2 = max(err2);
maxerr3 = max(err3);

fprintf('Precorrector: %d, %d, %d\n', maxerr1, maxerr2, maxerr3);

if((maxerr1>1e-2)|(maxerr2>1e-2)|(maxerr3>1e-2))
    disp('Error in Precorrecotr!')
end
delta1_0 = zeros(N,1);
delta2_0 = zeros(N,1);
delta1_0(2:N) = x1mat0(1:N-1,p+2)-x1vec0(2:N);
delta2_0(2:N) = x2mat0(1:N-1,p+2)-x2vec0(2:N);

    


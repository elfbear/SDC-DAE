function ploterror(tvec, errx1,errx2,erry)
figure(1)
hold on
plot(tvec(2:end),errx1,'r:+');
plot(tvec(2:end),errx2,'b:o');
plot(tvec(2:end),erry,'k:d');

xlabel('t') 
ylabel('Absolute error')
title('Error on coarse grids with 7 SDC iterations and p=5 for nonliear eqn')
legend('x1','x2','y');

hold off
print -dpdf n1.pdf
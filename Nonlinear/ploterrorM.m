function ploterrorM(mvec, errx1,errx2,erry)
figure(1)
hold on
plot(mvec,errx1,'r:+');
plot(mvec,errx2,'b:o');
plot(mvec,erry,'k:d');

xlabel('M') 
ylabel('-log10(|error|)')
title('Error on T with 3-15 SDC iterations and p=5 for nonliear case')
legend('x1','x2','y');

hold off
print -dpdf n2.pdf
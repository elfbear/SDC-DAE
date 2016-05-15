function ploterrorp(pvec, errx1,errx2,erry)
figure(1)
hold on
plot(pvec,errx1,'r:+');
plot(pvec,errx2,'b:o');
plot(pvec,erry, 'k:d');

xlabel('p') 
ylabel('-log10(|error|)')
title('Error on T with 5 SDC iterations and p=3-10 for nonliear case')
legend('x1','x2','y');

hold off
print -dpdf n3.pdf
nlsFILE = fopen('results8.txt','r');
fdmFILE = fopen('results8_fdm.txt','r');
nls = fscanf(nlsFILE, '%f %f %f', [3 Inf]);
fdm = fscanf(fdmFILE, '%f %f %f', [3 Inf]);

figure(1)
plot(nls(2,:),nls(3,:), '-o');
hold on
plot(fdm(2,:),fdm(3,:), '-o');
legend('NonLinear Shooting Method','NonLinear Finite Difference Method')
xlabel('x');
ylabel('Approximation')
title('Shooting Method and Finite Difference Method Solutions of BVP','FontSize',10)

figure(2)
plot(nls(2,:), nls(3,:)-fdm(3,:));
legend('Difference in Methods');
xlabel('x');
ylabel('Nonlinear Shooting - Finite Difference');
title('Nonlinear Shooting Method solution minus Finite Difference Method solution','FontSize',8);
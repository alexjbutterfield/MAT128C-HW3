%%% Problem 5.11.9

% set universal parameters
f = @(t,w) [32*w(1)+66*w(2)+2*t./3+2/3 -66*w(1)-133*w(2)-t./3-1/3];
inits = [1/3  1/3];
a = 0;
b = 0.5;

% set actual solution
sol = @(t) [2*t./3+2*exp(-t)./3-exp(-100*t)./3 -t./3-exp(-t)./3+2*exp(-100*t)./3];

% plot solution to figure 1
figure(1)
fplot(sol, [0 0.5]);
hold on

%% Part (a)
% accumulate data for part (a)
N = round((b-a)./0.1);
[t_a,w_a] = RK4SystemSolve(f,a,b,inits,N);
err_a = abs(sol(t_a) - w_a);

% plot data for part (a)
figure(2)
plot(t_a, log10(err_a));
legend('log_{10}(error_1)','log_{10}(error_2)');
xlabel('t_i');
ylabel('log_{10}(error(t_i))');
title('Plot of log_{10}(error) of RK4 solution with h = 0.1');

%% Part (b)
% accumulate data for part (b)
N = round((b-a)./0.025);
[t_b,w_b] = RK4SystemSolve(f,a,b,inits,N);
err_b = abs(sol(t_b) - w_b);

% add data to figure 1 to compare with actual solution
figure(1)
plot(t_b, w_b);
legend('u_1(t)','u_2(t)','w_{1}(t_i)','w_{2}(t_i)');
xlabel('t');
ylabel('u(t) and w(t)');
title('Comparision of actual solution and RK4 solution with h=0.025');
%%% Problem 11.1.2

% set universal parameters
p = @(x) 1;
q = @(x) 2;
r = @(x) cos(x);
a = 0;
b = pi./2;
alpha = -0.3;
beta = -0.1;

% set actual solutions with sol2 = d/dx (sol1) 
sol1 = @(x) -(sin(x) +3*cos(x))./10;
sol2 = @(x) -(cos(x) -3*sin(x))./10;

%% Part (a)

h = pi./4;
[ta,w1a,w2a] = LinearShootingMethod(p,q,r,a,b,alpha,beta,h);

figure(1)
plot(ta,[w1a w2a]);
hold on
fplot(sol1,[0 pi./2]);
fplot(sol2,[0 pi./2]);
legend('w_1(x_i)','w_2(x_i)','y(x)','dy(x)/dx');
xlabel('x_i');
ylabel('w_n(x_i) and y^{n}(x)');
title('Comparision of actual solution and LinearShooting solution with h= \pi/4');
hold off

TaW1 = table(ta',w1a,sol1(ta'), abs(w1a-sol1(ta')),'VariableNames',{'x','W1','y','err1'});
writetable(TaW1,'11-1-2aW1.csv');

TaW2 = table(ta',w2a,sol2(ta'),abs(w2a-sol2(ta')),'VariableNames',{'x','W2','Dy','err2'});
writetable(TaW2,'11-1-2aW2.csv');

%% Part (b)

h = pi./8;
[tb,w1b,w2b] = LinearShootingMethod(p,q,r,a,b,alpha,beta,h);

figure(2)
plot(tb,[w1b w2b]);
hold on
fplot(sol1,[0 pi./2]);
fplot(sol2,[0 pi./2]);
legend('w_1(x_i)','w_2(x_i)','y(x)','dy(x)/dx');
xlabel('x_i');
ylabel('w_n(x_i) and y^{n}(x)');
title('Comparision of actual solution and LinearShooting solution with h= \pi/8');
hold off

TbW1 = table(tb',w1b,sol1(tb'), abs(w1b-sol1(tb')),'VariableNames',{'x','W1','y','err1'});
writetable(TbW1,'11-1-2bW1.csv');

TbW2 = table(tb',w2b,sol2(tb'),abs(w2b-sol2(tb')),'VariableNames',{'x','W2','Dy','err2'});
writetable(TbW2,'11-1-2bW2.csv');
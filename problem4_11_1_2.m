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
hold off

%% Part (b)

h = pi./8;
[tb,w1b,w2b] = LinearShootingMethod(p,q,r,a,b,alpha,beta,h);

figure(2)
plot(tb,[w1b w2b]);
hold on
fplot(sol1,[0 pi./2]);
fplot(sol2,[0 pi./2]);
hold off

% FIXME: Add legends/axis-labels/title to figures and a few more comments to code
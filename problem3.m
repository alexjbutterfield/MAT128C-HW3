%%% COMPUTER PROJECT 1

% Set alpha
alpha = 0.01;

% Define animal system with added input for alpha
dadt = @(t,a) [2*a(1)-alpha*a(1).*a(2); -a(2)+alpha.*a(1).*a(2)];

%% Part (a)

% Has no MATLAB calculations required

%% Part (b)

% Generate initial a0's
test_a0 = [3000 50; 300 150; 150 1000; 102 198; 100 200];

i = 1;
for a0 = test_a0'
    [t,y] = ode45(dadt,[0 10],a0);
    
    figure(i)
    plot(t,y(:,1),'-o',t,y(:,2),'-o');
    legend('Rabbits','Foxes');
    xlabel('Time');
    ylabel('Population');
    title(['Lotka-Volterra System time solution for r_0 = ' num2str(a0(1)) ' and f_0 = ' num2str(a0(2))], 'FontSize', 10);
    
    figure(i+1)
    plot(y(:,1),y(:,2));
    xlabel('Rabbits');
    ylabel('Foxes');
    title(['Lotka-Volterra System phase space solution for r_0 = ' num2str(a0(1)) ' and f_0 = ' num2str(a0(2))], 'FontSize', 10);
    
    i = i+2;
end

%% Part (c)

R = 400;
dadt_mod = @(t,a) [2*(1-a(1)./R).*a(1)-0.01.*a(1).*a(2); -a(2)+0.01.*a(1).*a(2)];

a0 = test_a0(2,:);

[t,y] = ode45(dadt_mod,[0 10],a0);

figure(11)
plot(t,y(:,1),'-o',t,y(:,2),'-o');
legend('Rabbits','Foxes');
xlabel('Time');
ylabel('Population');
title(['Modified Lotka-Volterra System time solution for r_0 = ' num2str(a0(1)) ' and f_0 = ' num2str(a0(2))], 'FontSize', 10);
    
figure(12)
plot(y(:,1),y(:,2));
xlabel('Rabbits');
ylabel('Foxes');
title(['Lotka-Volterra System phase space solution for r_0 = ' num2str(a0(1)) ' and f_0 = ' num2str(a0(2))], 'FontSize', 10);

i = 1;
while i < 13
   fig = figure(i);
   gca.TitleFontSizeMultiplier = 1;
   set(fig,'PaperPositionMode','auto');
   print(['problem3-Figure-' num2str(i) '.png'],'-dpng','-r0')

   i = i + 1;
end
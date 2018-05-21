function [x,w1,w2] = LinearShootingMethod(p,q,r,a,b,alpha,beta, h)
    N = round((b-a)./h);
    
    fu = @(x,y) [y(2) p(x).*y(2)+q(x).*y(1)+r(x)];
    [tu, u] = RK4SystemSolve(fu,a,b,[alpha 0],N);
    
    fv = @(x,y) [y(2) p(x).*y(2)+q(x).*y(1)];
    [tv, v] = RK4SystemSolve(fv,a,b,[0 1],N);
    
    w = zeros(N+1,2);
    w(1,1) = alpha;
    w(1,2) = (beta - u(N+1,1))./v(N+1,1);
    
    w(2:N+1,:) = u(2:N+1,:) + w(1,2).*v(2:N+1,:);
    x = a:h:b;
    
    w1 = w(:,1);
    w2 = w(:,2);
end
function [t,w] = RK4SystemSolve(f,a,b,inits,N)
    m = size(inits,2);
    t = zeros(N+1,1);
    w = zeros(N+1,m);
    
    h = (b - a)./N;
    t(1) = a;
    w(1,:) = inits;
    for i = 1:N
        k1 = h.*f(t(i),w(i,:));
        k2 = h.*f(t(i) + h./2, w(i,:)+k1./2);
        k3 = h.*f(t(i) + h./2, w(i,:)+k2./2);
        k4 = h.*f(t(i)+h, w(i,:)+k3);
        w(i+1,:) = w(i,:) + (k1 + 2*k2 + 2*k3 + k4)./6;
        t(i+1) = t(i) + h;
    end
end
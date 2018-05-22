p = @(x) -2./x;
q = @(x) 0;
r = @(x) 0;
R1 = 2;
R2 = 4;
V1 = 110;
Vf = 0;

%% Part (a)

approx = zeros(4,1);

i = 1;
for h = [1 0.5 0.25 0.125]
    [t,w1,w2] = LinearShootingMethod(p,q,r,R1,R2,V1,Vf,h);
    approx(i) = w1((size(w1,1)+1)./2);
    i = i + 1;
end

%% Part (b)

solu = @(x) (V1.*R1./x).*(R2 - x)./(R2 - R1);

err = approx- solu(3);
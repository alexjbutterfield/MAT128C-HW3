p = @(x) -2./x;
q = @(x) 0;
r = @(x) 0;
R1 = 2;
R2 = 4;
V1 = 110;
Vf = 0;

%% Part (a)

h = 0.25;
[t,w1,w2] = LinearShootingMethod(p,q,r,R1,R2,V1,Vf,h);
w1(5)

%% Part (b)

solu = @(x) (V1.*R1./x).*(R2 - x)./(R2 - R1);
solu(3)
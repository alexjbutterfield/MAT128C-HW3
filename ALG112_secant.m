clear;
close all;

% Parameters
NAME = 'results4a_COMP2.txt';
A = 1; B = 2;
ALPHA = 1/2; BETA = 1/3;
N = 10;
NN = 100; TOL = 0.0001;
Y = @(x) (x+1).^(-1);
F = @(x,y,yp) y.^3 - y.*yp;
FY = @(x,y,yp) 3*y.^2 - yp;
FYP = @(x,y,yp) -y;
OUP = fopen(NAME,'wt');

% STEP 1
TK1 = (BETA-ALPHA)/(B-A);
W1 = zeros(1,N+1);
W2 = zeros(1,N+1);
H = (B-A)/N;
K = 1;
% TK already computed
OK = false;
% STEP 2
while K <= NN && OK == false 
    % STEP 3
    W1(1) = ALPHA;
    W2(1) = TK1;
    % STEP 4
    % Rung-Kutta method for systems is used in STEPS 5 and 6
    for I = 1 : N 
        %  STEP 5
        X = A+(I-1)*H;
        T = X+0.5*H;
        % STEP 6
        K11 = H*W2(I);
        K12 = H*F(X,W1(I),W2(I));
        K21 = H*(W2(I)+0.5*K12);
        K22 = H*F(T,W1(I)+0.5*K11,W2(I)+0.5*K12);
        K31 = H*(W2(I)+0.5*K22);
        K32 = H*F(T,W1(I)+0.5*K21,W2(I)+0.5*K22);
        K41 = H*(W2(I)+K32);
        K42 = H*F(X+H,W1(I)+K31,W2(I)+K32);
        W1(I+1) = W1(I)+(K11+2*(K21+K31)+K41)/6;
        W2(I+1) = W2(I)+(K12+2*(K22+K32)+K42)/6;
    end
    
    
    % STEP 7
    % test for accuracy
    if abs(W1(N+1)-BETA) < TOL 
        % STEP 8
        I = 0;
        xvalues = zeros(N,1);
        yvalues = zeros(N,1);
        ytrues = zeros(N,1);
        fprintf(OUP, '%3s %13s %13s %13s %13s %13s\n', 'N', 'X', 'W1', 'Y', 'Error', 'W2');
        fprintf(OUP, '%3d %13.8f %13.8f %13.8f %13.8f %13.8f\n', I, A, ALPHA, Y(A), abs(ALPHA - Y(A)), TK1);
        for I = 1 : N 
            J = I+1;
            X = A+I*H;
            xvalues(I) = X;
            yvalues(I) = W1(J);
            ytrues(I) = Y(X);
            fprintf(OUP, '%3d %13.8f %13.8f %13.8f %13.8f %13.8f\n', I, X, W1(J), Y(X), abs(Y(X)-W1(J)), W2(J));
        end
        fprintf(OUP, 'Convergence in %d iterations\n', K);
        fprintf(OUP, ' t = %14.7e\n', TK1);
        % STEP 9
        OK = true;
    else
        % STEP 10
        % Newton's method applied to improve TK
        
        if K == 1
            TK0 = TK1;
            TK1 = TK1 + (BETA - W1(N+1))./(B - A);
            W1_TK0 = W1(N+1);
        else
            TKbuff = TK1;
            TK1 = TK1 - (W1(N+1) - BETA).*(TK1 - TK0)./(W1(N+1) - W1_TK0);
            TK0 = TK1;
            W1_TK0 = W1(N+1);
        end
        K = K+1;
    end
end
% STEP 11

% method failed
if OK == false 
    fprintf(OUP, 'Method failed after %d iterations\n', NN);
end


if OUP ~= 1 
    fclose(OUP);
    fprintf(1,'Output file %s created successfully \n',NAME);
end

figure;
plot(xvalues, yvalues, 'ro'); hold on;
plot(xvalues, ytrues, 'b-');
legend('W1', 'Y')

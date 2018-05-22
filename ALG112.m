% NONLINEAR SHOOTING ALGORITHM 11.2
%
% To approximate the solution of the nonlinear boundary-value problem
%
%          Y'' = F(X,Y,Y'), A<=X<=B, Y(A) = ALPHA, Y(B) = BETA:
%
%
% INPUT:   Endpoints A,B; boundary conditions ALPHA, BETA; number of
%          subintervals N; tolerance TOL; maximum number of iterations M.
%
% OUTPUT:  Approximations W(1,I) TO Y(X(I)); W(2,I) TO Y'(X(I))
%          for each I=0,1,...,N or a message that the maximum
%          number of iterations was exceeded.
clear;
close all;

% Parameters
NAME = 'results2.txt';
A = -1; B = 0;
ALPHA = 1/2; BETA = 1/3;
N = 4;
NN = 100; TOL = 0.0001;
Y = @(x) 1/(x+3);
F = @(x,y,yp) 2*y^3;
FY = @(x,y,yp) 6*y^2;
FYP = @(x,y,yp) 0;
OUP = fopen(NAME,'wt');

% STEP 1
TK = (BETA-ALPHA)/(B-A);
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
    W2(1) = TK;
    U1 = 0 ;
    U2 = 1;
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
        K11 = H*U2;
        K12 = H*(FY(X,W1(I),W2(I))*U1+FYP(X,W1(I),W2(I))*U2);
        K21 = H*(U2+0.5*K12);
        K22 = H*(FY(T,W1(I),W2(I))*(U1+0.5*K11)+FYP(T,W1(I),W2(I))*(U2+0.5*K21));
        K31 = H*(U2+0.5*K22);
        K32 = H*(FY(T,W1(I),W2(I))*(U1+0.5*K21)+FYP(T,W1(I),W2(I))*(U2+0.5*K22));
        K41 = H*(U2+K32);
        K42 = H*(FY(X+H,W1(I),W2(I))*(U1+K31)+FYP(X+H,W1(I),W2(I))*(U2+K32));
        U1 = U1+(K11+2*(K21+K31)+K41)/6;
        U2 = U2+(K12+2*(K22+K32)+K42)/6;
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
        fprintf(OUP, '%3d %13.8f %13.8f %13.8f %13.8f %13.8f\n', I, A, ALPHA, Y(A), abs(ALPHA - Y(A)), TK);
        for I = 1 : N 
            J = I+1;
            X = A+I*H;
            xvalues(I) = X;
            yvalues(I) = W1(J);
            ytrues(I) = Y(X);
            fprintf(OUP, '%3d %13.8f %13.8f %13.8f %13.8f %13.8f\n', I, X, W1(J), Y(X), abs(Y(X)-W1(J)), W2(J));
        end
        fprintf(OUP, 'Convergence in %d iterations\n', K);
        fprintf(OUP, ' t = %14.7e\n', TK);
        % STEP 9
        OK = true;
    else
        % STEP 10
        % Newton's method applied to improve TK
        TK = TK-(W1(N+1)-BETA)/U1;
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

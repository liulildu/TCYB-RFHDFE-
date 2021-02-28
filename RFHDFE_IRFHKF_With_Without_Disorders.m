% A distributed fusion estimation approach based on the robust finite horizon filtering for estimated results with or without packet disorders.
clear;      clc;            close all;
T=0.1;    A=[0.9  T  T^2/2;  0  0.9  T;  0  0  0.9];
B=[T^2/2; T; 1];        C=[0.6  0.8  1];        N=5;        % Largest delays
F1=[0.1; 0.1; 0.1];     H1=0.8;     E=[0.02  0.02  0.02];        iter=305;
beta=2; a=3;    F=zeros(iter,1);
Peta=0.09;      Q_k=Peta;       R_k=beta*Q_k*beta';     S_k=Q_k*beta';
w=sqrt(Q_k)*randn(iter,1);      v=sqrt(R_k)*randn(iter,1);  % Actual State
w1=sqrt(Q_k)*randn(iter,1);     v1=sqrt(R_k)*randn(iter,1); % Delay State
% x(k|k-1);                 % x(k|k);               % x(k|t1);
x1=zeros(3,iter);           x2=zeros(3,iter);       x3=zeros(3,iter);             
P1=zeros(3,3*iter);           % P1=x(k)*x(k)';
Trtheta1=zeros(iter,1);             % Trace(theta1);
% sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))') % theta1=E((x(k)-x(k|k))(x(k)-x(k|k))')
sigma1=zeros(3,3*iter);         theta1=zeros(3,3*iter);
% lambda1=sigma1*C'           delta1=A*sigma1           Xi1=C*sigma1...
lambda1=zeros(3,iter);        delta1=zeros(3,iter);     Xi1=zeros(iter,1);
% M1=1/a-E*P*E'                     M2=1/a-E*sigma1*E'
M1=zeros(iter,1);                   M2=zeros(iter,1);
tru=zeros(3,iter);                  z=zeros(iter,1);    % z=(C+HFE)tru(k)+v_k 
tru1=zeros(3,iter);                 z1=zeros(iter,1);   % z1=(C+HFE)tru(k)+v_k
C1=zeros(1,3*iter);     K1=zeros(3,iter);    
A1=zeros(3,3*iter);     L1=zeros(3,iter);       % Filtering parameter
tau=zeros(iter,1);      t=zeros(iter,1);        % Transmission delay
tau1=zeros(iter,1);     t1=zeros(iter,1);       % Received packet with time-stamped
s1=zeros(iter,1);       % Re-organized time delay
b1=zeros(3,iter);       b2=zeros(3,iter);     b3=zeros(3,iter);        % eig

P=0.01*eye(3); 
tru(:,1)=[1;  1;  1];            % ¼ÇÂ¼ÕæÊµ×´Ì¬
x1(:,1)=tru(:,1);               % x(k|k-1);
x2(:,1)=tru(:,1);               % x(k|k);
sigma1(:,1:3)=P;                % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))')
P1(:,1:3)=tru(:,1)*tru(:,1)'+P; % P1=x(k)*x(k);
% b2(:,1)=eig(sigma1(:,1:3));     b3(:,1)=eig(P1(:,1:3));

for k=1:iter
    if k<=N
        tau(k,1)=round(rand(1,1)*k);
    else
        tau(k,1)=round(rand(1,1)*N);
    end
    F(k,1)=sin(0.6*k); 
    M1_IS(k,1)=(1/a)-E*P1(:,3*k-2:3*k)*E';       % M1=1/a-E*P*E'
    M2_IS(k,1)=(1/a)-E*sigma1(:,3*k-2:3*k)*E';   % M2=1/a-E*sigma1*E'
    z(k,1)=(C+H1*F(k,1)*E)*tru(:,k)+v(k,1);
    tru(:,k+1)=(A+F1*F(k,1)*E)*tru(:,k)+B*w(k,1);
end

[sigma1, Trtheta1, P1, x3, M1, M2, fv, et_ZOH]=RFHDFE_ZOH_Function_TCYB(T, A, B, C, E, F1, H1, F, a, beta,...
    w1, v1, Q_k, R_k, S_k, tru, z, x1, x2, sigma1, P1, iter, tau, N);
[sigmaL1, TrthetaL1, PL1, xL3, ML1, ML2, fvL, et_LZOH]=RFHDFE_LZOH_Function_TCYB(T, A, B, C, E, F1, H1, F, a, beta,...
    w1, v1, Q_k, R_k, S_k, tru, z, x1, x2, sigma1, P1, iter, tau, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  RFHKF_Function  with time-delay   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sigma11, Trtheta11, P11, x22, MIS1, MIS2, fvIS, et_IS]=IRFHKF_Function_Delay_TCYB(T, A, B, C, E, F1, H1, F, a, beta,...
    Q_k, R_k, S_k, M1_IS, M2_IS, tru,  z, x1,sigma1, P1, iter, tau, N);

d1=zeros(3,3);
for k=1:iter
    b1(1,k)=fv(1,3*k-2);
    b1(2,k)=fv(2,3*k-1);
    b1(3,k)=fv(3,3*k);
    b2(1,k)=fvL(1,3*k-2);
    b2(2,k)=fvL(2,3*k-1);
    b2(3,k)=fvL(3,3*k);
    b3(1,k)=fvIS(1,3*k-2);
    b3(2,k)=fvIS(2,3*k-1);
    b3(3,k)=fvIS(3,3*k);
    d1(1,1)=d1(1,1)+b1(1,k);
    d1(1,2)=d1(1,2)+b1(2,k);
    d1(1,3)=d1(1,3)+b1(3,k);
    d1(2,1)=d1(2,1)+b2(1,k);
    d1(2,2)=d1(2,2)+b2(2,k);
    d1(2,3)=d1(2,3)+b2(3,k);
    d1(3,1)=d1(3,1)+b3(1,k);
    d1(3,2)=d1(3,2)+b3(2,k);
    d1(3,3)=d1(3,3)+b3(3,k);
end
d1(1,1)=d1(1,1)/iter;
d1(1,2)=d1(1,2)/iter;
d1(1,3)=d1(1,3)/iter;
d1(2,1)=d1(2,1)/iter;
d1(2,2)=d1(2,2)/iter;
d1(2,3)=d1(2,3)/iter;
d1(3,1)=d1(3,1)/iter;
d1(3,2)=d1(3,2)/iter;
d1(3,3)=d1(3,3)/iter;

iter=300;

%%%%%%%%%%%%%%% Position states %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;             % Fig 1
plot(tru(1,1:iter),'r');        hold on;    
plot(x3(1,1:iter),'c');         hold on; 
plot(xL3(1,1:iter),'b');        hold on;
plot(x22(1,1:iter),'g');        hold on;
xlabel('t/step');     ylabel('States');
legend('Actual state of position','Estimated state with packet disorders',...
    'Estimated state of RFHDEF','Estimated state of IRFHKF');
%%%%%%%%%%%%%% Velocity states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;             % Fig 2
plot(tru(2,1:iter),'r');        hold on;    
plot(x3(2,1:iter),'c');         hold on;
plot(xL3(2,1:iter),'b');        hold on;
plot(x22(2,1:iter),'g');        hold on;
xlabel('t/step');     ylabel('States');
legend('Actual state of velocity','Estimated state with packet disorders',...
    'Estimated state of RFHDFE','Estimated state of IRFHKF');
%%%%%%%%%%%%%%% Acceleration states
figure;             % Fig 3
plot(tru(3,1:iter),'r');        hold on;    
plot(x3(3,1:iter),'c');         hold on;
plot(xL3(3,1:iter),'b');        hold on;
plot(x22(3,1:iter),'g');        hold on;
xlabel('t/step');     ylabel('States');
legend('Actual state of acceleration','Estimated state with packet disorders',...
    'Estimated state of RFHDFE', 'Estimated state of IRFHKF');

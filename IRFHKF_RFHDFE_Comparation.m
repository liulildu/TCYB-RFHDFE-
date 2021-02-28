% A distributed fusion estimation approach based on the robust finite horizon filtering
clear;      clc;            close all;
T=0.1;    A=[0.9  T  T^2/2;  0  0.9  T;  0  0  0.9]; L=3;   N=5;  % Largest delays
B=[T^2/2; T; 1];        C=[0.6,  0.8,  1; 1,  0.8,  0.5; 0.3,  1, 0.7];       
F1=[0.1; 0.1; 0.1];     H1=0.8;     E=[0.02  0.02  0.02];        iter=300;
beta=[2; 0.8;  1]; a=3;    F=zeros(iter,1);
Peta=0.09;      Q_k=Peta;       R_k=beta*Q_k*beta';     S_k=beta*Q_k;
w=sqrt(Q_k)*randn(iter,1);      
v(:,1)=randn(iter,1)*sqrt(R_k(1,1));  v(:,2)=randn(iter,1)*sqrt(R_k(2,2));  v(:,3)=randn(iter,1)*sqrt(R_k(3,3));% Actual State
w1=sqrt(Q_k)*randn(iter,1);     
v1(:,1)=randn(iter,1)*sqrt(R_k(1,1));  v1(:,2)=randn(iter,1)*sqrt(R_k(2,2));  v1(:,3)=randn(iter,1)*sqrt(R_k(3,3));% Delay State
% x(k|k-1);                 % x(k|k);               % x(k|t1);
x1=zeros(3*L,iter);         x2=zeros(3*L,iter);     x3=zeros(3*L,iter);             
P1=zeros(3,3*iter);           % P1=x(k)*x(k)';
Trtheta1=zeros(iter,L);             % Trace(theta1);
% sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))') % theta1=E((x(k)-x(k|k))(x(k)-x(k|k))')
sigma1=zeros(3*L,3*iter);         theta1=zeros(3*L,3*iter);
% lambda1=sigma1*C'           delta1=A*sigma1           Xi1=C*sigma1...
lambda1=zeros(3*L,iter);      delta1=zeros(3*L,iter);   Xi1=zeros(iter,1*L);
% M1=1/a-E*P*E'                     M2=1/a-E*sigma1*E'
M1=zeros(iter,1);                   M2=zeros(iter,1*L);
tru=zeros(3,iter);                  z=zeros(iter,1*L);    % z=(C+HFE)tru(k)+v_k 
tru1=zeros(3*L,iter);               z1=zeros(iter,1*L);   % z1=(C+HFE)tru(k)+v_k
C1=zeros(1*L,3*iter);     K1=zeros(3*L,iter);    
A1=zeros(3*L,3*iter);     L1=zeros(3*L,iter);       % Filtering parameter
tau=zeros(iter,1);      t=zeros(iter,1);        % Transmission delay
tau1=zeros(iter,1);     t1=zeros(iter,1);       % Received packet with time-stamped
s1=zeros(iter,1);       % Re-organized time delay
b1=zeros(3,iter);       b2=zeros(3,iter);       b3=zeros(3,iter);   % eig
P=0.01*eye(3); 
tru(:,1)=[1;  1;  1];            % actual state
x1(1:3,1)=tru(:,1);  x1(4:6,1)=tru(:,1);   x1(7:9,1)=tru(:,1);             % x(k|k-1);
for i=1:L
    for j=1:L
        sigma1(3*i-2:3*i,3*j-2:3*j)=P;  % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))')
    end
end 
P1(:,1:3)=tru(:,1)*tru(:,1)'+P; % P1=x(k)*x(k);
%%%%%%%%% Time-delay, State and measurement
for k=1:iter
    M1(k,1)=(1/a)-E*P1(:,3*k-2:3*k)*E';       % M1=1/a-E*P*E'
    F(k,1)=sin(0.6*k);
    tru(:,k+1)=(A+F1*F(k,1)*E)*tru(:,k)+B*w(k,1);
    for i=1:L
        M2(k,i)=(1/a)-E*sigma1(3*i-2:3*i,3*k-2:3*k)*E';   % M2=1/a-E*sigma1*E'
        z(k,i)=(C(i,:)+H1*F(k,1)*E)*tru(:,k)+v(k,i);
    end
end

[sigma1, Trtheta1, P1, x2, Trtheta2]=RFHDFE_3D_Distributed_Function(T, A, B, C, E, F1, H1, F, a, beta,...
    Q_k, R_k, S_k, M1, M2, tru,  z, x1,sigma1, P1, iter, L);
[sigma11, Trtheta11, P11, x22]=IRFHKF_Function(T, A, B, C(1,:), E, F1, H1, F, a, beta(1,1),...
    Q_k, R_k(1,1), S_k(1,1), M1, M2(:,1), tru,  z(:,1), x1(1:3,:),sigma1(1:3,:), P1, iter);

figure;
plot(Trtheta1(:,1),'Color',[0.5 0.5 0.5]);   hold on;
plot(Trtheta1(:,2),'c');   hold on;
plot(Trtheta1(:,3),'Color',[1 0.5 0]);   hold on;
plot(Trtheta2(:,1),'r');   hold on;
plot(Trtheta11(:,1),'b');   hold on;
xlabel('t/step');     ylabel('Trace');
legend('Trace of error covariance for sensor 1','Trace of error covariance for sensor 2',...
    'Trace of error covariance for sensor 3','Trace of error covariance of fusion estimation','Trace of IRFHKF');
text('Interpreter','latex','String','$\zeta_1  = 2$');
text('Interpreter','latex','String','$\zeta_2  = 0.8$');
text('Interpreter','latex','String','$\zeta_3  = 1$');

for k=1:iter    
    c1(1,k)=sigma1(1, 3*k-2);
    c1(2,k)=sigma1(2,3*k-1);    
    c1(3,k)=sigma1(3,3*k);    
    c1(4,k)=sigma1(4,3*k-2);  
    c1(5,k)=sigma1(5,3*k-1);    
    c1(6,k)=sigma1(6,3*k);     
    c1(7,k)=sigma1(7,3*k-2);   
    c1(8,k)=sigma1(8,3*k-1);    
    c1(9,k)=sigma1(9,3*k);
   
    c2(1,k)=sigma11(1,3*k-2);
    c2(2,k)=sigma11(2,3*k-1);
    c2(3,k)=sigma11(3,3*k);
end
d1(1,1)=min(c1(1,:));   d1(1,2)=max(c1(1,:));
d1(1,3)=min(c1(2,:));   d1(1,4)=max(c1(2,:));
d1(1,5)=min(c1(3,:));   d1(1,6)=max(c1(3,:));
d1(1,7)=min(Trtheta1(:,1));   d1(1,8)=max(Trtheta1(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1(2,1)=min(c1(4,:));   d1(2,2)=max(c1(4,:));
d1(2,3)=min(c1(5,:));   d1(2,4)=max(c1(5,:));
d1(2,5)=min(c1(6,:));   d1(2,6)=max(c1(6,:));
d1(2,7)=min(Trtheta1(:,2));   d1(2,8)=max(Trtheta1(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1(3,1)=min(c1(7,:));   d1(3,2)=max(c1(7,:));
d1(3,3)=min(c1(8,:));   d1(3,4)=max(c1(8,:));
d1(3,5)=min(c1(9,:));   d1(3,6)=max(c1(9,:));
d1(3,7)=min(Trtheta1(:,3));   d1(3,8)=max(Trtheta1(:,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1(4,1)=min(c2(1,:));   d1(4,2)=max(c2(1,:));
d1(4,3)=min(c2(2,:));   d1(4,4)=max(c2(2,:));
d1(4,5)=min(c2(3,:));   d1(4,6)=max(c2(3,:));
d1(4,7)=min(Trtheta11(:,1));   d1(4,8)=max(Trtheta11(:,1));

%%%%%%%%%%%%%%% Position states %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(tru(1,1:iter),'r');       hold on;    
plot(x2(1,1:iter),'b');        hold on;
plot(x2(4,1:iter),'Color',[1 0.5 0]);        hold on;
plot(x2(7,1:iter),'Color',[0.5 0.5 0.5]);        hold on;
plot(x22(1,1:iter),'c');       hold on;

xlabel('t/step');     ylabel('States');
legend('Actual state of position','Estimated state of RFHDFE by sensor 1',...
    'Estimated state of RFHDFE by sensor 2','Estimated state of RFHDFE by sensor 3','Estimated state of IRFHKF');

%%%%%%%%%%%%%% Velocity states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(tru(2,1:iter),'r');       hold on;    
plot(x2(2,1:iter),'b');        hold on; 
plot(x2(5,1:iter),'Color',[1 0.5 0]);        hold on;
plot(x2(8,1:iter),'Color',[0.5 0.5 0.5]);        hold on;
plot(x22(2,1:iter),'c');       hold on;
xlabel('t/step');     ylabel('States');
legend('Actual state of velocity','Estimated state of RFHDFE by sensor 1',...
    'Estimated state of RFHDFE by sensor 2','Estimated state of RFHDFE by sensor 3','Estimated state of IRFHKF');

%%%%%%%%%%%%%%% Acceleration states
figure;
plot(tru(3,1:iter),'r');       hold on;    
plot(x2(3,1:iter),'b');        hold on;
plot(x2(6,1:iter),'Color',[1 0.5 0]);        hold on;
plot(x2(9,1:iter),'Color',[0.5 0.5 0.5]);        hold on;
plot(x22(3,1:iter),'c');       hold on;
xlabel('t/step');     ylabel('States');
legend('Actual state of acceleration','Estimated state of RFHDFE by sensor 1',...
    'Estimated state of RFHDFE by sensor 2','Estimated state of RFHDFE by sensor 3','Estimated state of IRFHKF');

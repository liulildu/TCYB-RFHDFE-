% A distributed fusion estimation approach based on the robust finite horizon filtering 
clear;      clc;            close all;
T=0.1;    A=[0.9  T  T^2/2;  0  0.9  T;  0  0  0.9]; L=3;   N=5;  % Largest delays
B=[T^2/2; T; 1];        C=[0.6,  0.8,  1; 1,  0.8,  0.5; 0.3,  1, 0.7];       
F1=[0.1; 0.1; 0.1];     H1=0.8;     E=[0.02  0.02  0.02];        iter=305;
beta=[2; 0.8;  1]; a=3;    F=zeros(iter,1);     st=[0  0; 0  0; 0   0; 0   0]; %First: cuptime; Second: end time 
Peta=0.09;      Q_k=Peta;       R_k=beta*Q_k*beta';     S_k=beta*Q_k;
w=sqrt(Q_k)*randn(iter,1);              v(:,1)=randn(iter,1)*sqrt(R_k(1,1));  
v(:,2)=randn(iter,1)*sqrt(R_k(2,2));    v(:,3)=randn(iter,1)*sqrt(R_k(3,3));% Actual State
w1=sqrt(Q_k)*randn(iter,1);             v1(:,1)=randn(iter,1)*sqrt(R_k(1,1));
v1(:,2)=randn(iter,1)*sqrt(R_k(2,2));   v1(:,3)=randn(iter,1)*sqrt(R_k(3,3));% Delay State
% x(k|k-1);                 % x(k|k);               % x(k|t1);
x1=zeros(3*L,iter);         x2=zeros(3*L,iter);     x3=zeros(3*L,iter);     x4=zeros(3,iter);            
P1=zeros(3,3*iter);             % P1=x(k)*x(k)';
% Trace(theta1);                % Trace(Pi)
Trtheta1=zeros(iter,L);         Trtheta2=zeros(iter,1);          
% sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))') % theta1=E((x(k)-x(k|k))(x(k)-x(k|k))')
sigma1=zeros(3*L,3*iter);                   theta1=zeros(3*L,3*iter);
% lambda1=sigma1*C'           delta1=A*sigma1           Xi1=C*sigma1...
lambda1=zeros(3*L,iter);      delta1=zeros(3*L,iter);   Xi1=zeros(iter,1*L);
% sigma1=E((x(k)-x_i(k|k-1))(x(k)-x_j(k|k-1))')     % Pi1=E((x(k)-x(k|t))(x(k)-x(k|t))')
sigma2=zeros(3*L,3*L*iter);                           Pi1=zeros(3,3*iter);
Omega1=zeros(3,3*L*iter);       fv=zeros(3*L,3*iter);     %fv=(tru(:,k)-x3(3*i-2:3*i,k))*(tru(:,k)-x3(3*i-2:3*i,k))';
b1=zeros(4*L,iter);               d1=zeros(4*L,1);
% M1=1/a-E*P*E'                 M2=1/a-E*sigma1*E' 
M1=zeros(iter,1);               M2=zeros(iter,1*L);     M2(1,1:3)=[1  1  1];
tru=zeros(3,iter);              z=zeros(iter,1*L);    % z=(C+HFE)tru(k)+v_k 
tru1=zeros(3*L,iter);           z1=zeros(iter,1*L);   % z1=(C+HFE)tru(k)+v_k
C1=zeros(1*L,3*iter);           K1=zeros(3*L,iter);    
A1=zeros(3*L,3*iter);           L1=zeros(3*L,iter);       % Filtering parameter
tau=zeros(iter,1);      t=zeros(iter,1);        % Transmission delay
tau1=zeros(iter,1);     t1=zeros(iter,1);       % Received packet with time-stamped
s1=zeros(iter,1);       % Re-organized time delay  
P=0.01*eye(3);          tru(:,1)=[1;  1;  1];            % actual state
TrP1=zeros(3,iter);             %Trace()
for i=1:L
    for j=1:L
        sigma1(3*i-2:3*i,3*j-2:3*j)=P;
    end
end
x1(1:3,1)=tru(:,1);     x1(4:6,1)=tru(:,1);     x1(7:9,1)=tru(:,1); % x(k|k-1);
%sigma1(1:3,1:3)=P;      sigma1(4:6,1:3)=P;      sigma1(7:9,1:3)=P;  % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))')
P1(:,1:3)=tru(:,1)*tru(:,1)'+P; % P1=x(k)*x(k);
I1=[eye(3)  eye(3)  eye(3)]';
%%%%%%%%% Time-delay, State and measurement
for k=1:iter
    if k<=N
        tau(k,1)=round(rand(1,1)*k);
    else
        tau(k,1)=round(rand(1,1)*N);
    end
    t(k,1)=k-tau(k,1);
    t1(k,1)=t(k,1);
    tau1(k,1)=k-t1(k,1);
    if k>1 && t1(k,1)<t1(k-1,1) 
        t1(k,1)=t1(k-1,1);
        tau1(k,1)=k-t1(k,1);
    end
    if k>1
        s1(k,1)=t1(k,1)-t1(k-1,1);
    end
    F(k,1)=sin(0.6*k);
    z(k,1)=(C(1,:)+H1*F(k,1)*E)*tru(:,k)+v(k,1);
    z(k,2)=(C(2,:)+H1*F(k,1)*E)*tru(:,k)+v(k,2);
    z(k,3)=(C(3,:)+H1*F(k,1)*E)*tru(:,k)+v(k,3);
    tru(:,k+1)=(A+F1*F(k,1)*E)*tru(:,k)+B*w(k,1);
end

for m=1:L
    st(m,1)=cputime; 
    for k=1:iter
        if t1(k,1)==0
            x2(3*m-2:3*m,t1(k,1)+1)=tru(:,1);
            x3(3*m-2:3*m,k)=tru(:,1);
            if tau1(k,1)>=1
                x3(3*m-2:3*m,k)=((N-tau1(k,1)+1)/N)*tru(:,1);
            end
            continue;
        else
            tru1(3*m-2:3*m,t1(k,1))=tru(:,t1(k,1));
            z1(t1(k,1),m)=z(t1(k,1),m);   % z1=(C+HFE)tru(k)+v_k
        end      
        if s1(k,1)<=1
            if s1(k,1)==1
                M1(t1(k,1),1)=(1/a)-E*P1(:,3*t1(k,1)-2:3*t1(k,1))*E';       % M1=1/a-E*P*E'
                M2(t1(k,1),m)=(1/a)-E*sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*E';   % M2=1/a-E*sigma1*E'
                lambda1(3*m-2:3*m,t1(k,1))=(eye(3)+sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M2(t1(k,1),m))*E)*...
                    sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*C(m,:)';% lambda1=sigma1*C'
                delta1(3*m-2:3*m,t1(k,1))=A*sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*...
                    (eye(3)+E'*inv(M2(t1(k,1),m))*E*sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1)))*C(m,:)'+...
                    (1/a)*F1*H1'+B*S_k(m,1);% delta1=A*sigma1
                % Xi1=C*sigma1
                Xi1(t1(k,1),m)=C(m,:)*sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*...
                    (eye(3)+E'*inv(M2(t1(k,1),m))*E*sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1)))*C(m,:)'+...
                    (1/a)*H1*H1'+R_k(m,m);
                % theta1=E((x(k)-x(k|k))(x(k)-x(k|k))')
                theta1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))=sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))+...
                    sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M1(t1(k,1),1))*E*sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))-...
                    lambda1(3*m-2:3*m,t1(k,1))*inv(Xi1(t1(k,1),1))*lambda1(3*m-2:3*m,t1(k,1))';
                Trtheta1(t1(k,1),m)=trace(theta1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1)));             % Trace(theta1);
                % Filtering parameters C1, K1, A1, L1 
                C1(m,3*t1(k,1)-2:3*t1(k,1))=C(m,:)*(eye(3)+...
                    sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M2(t1(k,1),m))*E);
                K1(3*m-2:3*m,t1(k,1))=lambda1(3*m-2:3*m,t1(k,1))*inv(Xi1(t1(k,1),m));    
                A1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))=A*(eye(3)+...
                    sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M2(t1(k,1),m))*E);
                L1(3*m-2:3*m,t1(k,1))=delta1(3*m-2:3*m,t1(k,1))*inv(Xi1(t1(k,1),m));
                x2(3*m-2:3*m,t1(k,1))=x1(3*m-2:3*m,t1(k,1))+K1(3*m-2:3*m,t1(k,1))*...
                    (z1(t1(k,1),m)-C1(m,3*t1(k,1)-2:3*t1(k,1))*x1(3*m-2:3*m,t1(k,1)));         % x2=x(k|k);
                x1(3*m-2:3*m,t1(k,1)+1)=A1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*x1(3*m-2:3*m,t1(k,1))+L1(3*m-2:3*m,t1(k,1))*...
                    (z1(t1(k,1),m)-C1(m,3*t1(k,1)-2:3*t1(k,1))*x1(3*m-2:3*m,t1(k,1)));           % x1=x(k|k-1);
                % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))') 
                sigma1(3*m-2:3*m,3*(t1(k,1)+1)-2:3*(t1(k,1)+1))=A*sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*...
                    (eye(3)+E'*inv(M2(t1(k,1),m))*E*sigma1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1)))*A'-...
                    delta1(3*m-2:3*m,t1(k,1))*inv(Xi1(t1(k,1),m))*delta1(3*m-2:3*m,t1(k,1))'+...
                    B*Q_k*B'+(1/a)*F1*F1';  
                % P1=x(k)*x(k)'; 
                P1(:,3*(t1(k,1)+1)-2:3*(t1(k,1)+1))=A*P1(:,3*t1(k,1)-2:3*t1(k,1))*A'+...
                    A*P1(:,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M1(t1(k,1),1))*E*P1(:,3*t1(k,1)-2:3*t1(k,1))*A'+...
                    (1/a)*F1*F1'+B*Q_k*B';
            end
            if tau1(k,1)==0
                x3(3*m-2:3*m,k)=x2(3*m-2:3*m,t1(k,1));
            elseif tau1(k,1)==1
                x3(3*m-2:3*m,k)=x1(3*m-2:3*m,t1(k,1)+1);
            else
                x3(3*m-2:3*m,k)=((N-tau1(k,1)+1)/N)*x1(3*m-2:3*m,t1(k,1)+1);
            end
        end
        % if s1(k,1)>1 || tau1(k,1)>1
        if s1(k,1)>1
            for i=1:s1(k,1)
                M1(t1(k-1,1)+i,1)=(1/a)-E*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E';       % M1=1/a-E*P*E'
                M2(t1(k-1,1)+i,m)=(1/a)-E*sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E';   % M2=1/a-E*sigma1*E'
                lambda1(3*m-2:3*m,t1(k-1,1)+i)=(eye(3)+sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*...
                    inv(M2(t1(k-1,1)+i,m))*E)*sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*C(m,:)';% lambda1=sigma1*C'
                delta1(3*m-2:3*m,t1(k-1,1)+i)=A*sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*...
                    (eye(3)+E'*inv(M2(t1(k-1,1)+i,m))*E*...
                    sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*C(m,:)'+(1/a)*F1*H1'+B*S_k(m,1);% delta1=A*sigma1
                % Xi1=C*sigma1
                Xi1(t1(k-1,1)+i,m)=C(m,:)*sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*...
                    (eye(3)+E'*inv(M2(t1(k-1,1)+i,m))*E*...
                    sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*C(m,:)'+(1/a)*H1*H1'+R_k(m,m);
                % theta1=E((x(k)-x(k|k))(x(k)-x(k|k))')
                theta1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))=sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))+...
                    sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(M1(t1(k-1,1)+i,1))*E*...
                    sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))-...
                    lambda1(3*m-2:3*m,t1(k-1,1)+i)*inv(Xi1(t1(k-1,1)+i,m))*lambda1(3*m-2:3*m,t1(k-1,1)+i)';
                Trtheta1(t1(k-1,1)+i,m)=trace(theta1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)));% Trace(P1);
                % Filtering parameters C1, K1, A1, L1 
                C1(m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))=C(m,:)*(eye(3)+...
                    sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(M2(t1(k-1,1)+i,m))*E);
                K1(3*m-2:3*m,t1(k-1,1)+i)=lambda1(3*m-2:3*m,t1(k-1,1)+i)*inv(Xi1(t1(k-1,1)+i,m));    
                A1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))=A*(eye(3)+...
                    sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(M2(t1(k-1,1)+i,m))*E);
                L1(3*m-2:3*m,t1(k-1,1)+i)=delta1(3*m-2:3*m,t1(k-1,1)+i)*inv(Xi1(t1(k-1,1)+i,m));
                if i<s1(k,1)
                    if t1(k-1,1)==0
                        t1(k-1,1)=t1(k-1,1)+1;
                        tru1(3*m-2:3*m,t1(k-1,1))=tru(:,1);
                        z1(t1(k-1,1),m)=z(1,m);
                        x2(3*m-2:3*m,t1(k-1,1))=x1(3*m-2:3*m,t1(k-1,1))+K1(3*m-2:3*m,t1(k-1,1))*...
                            (z1(t1(k-1,1),m)-C1(m,3*t1(k-1,1)-2:3*t1(k-1,1))*x1(3*m-2:3*m,t1(k-1,1)));         % x2=x(k|k);
                        x1(3*m-2:3*m,t1(k-1,1)+1)=A1(3*m-2:3*m,3*t1(k-1,1)-2:3*t1(k-1,1))*x1(3*m-2:3*m,t1(k-1,1))+...
                            L1(3*m-2:3*m,t1(k-1,1))*(z1(t1(k-1,1),m)-...
                            C1(m,3*t1(k-1,1)-2:3*t1(k-1,1))*x1(3*m-2:3*m,t1(k-1,1)));           % x1=x(k|k-1);
                        sigma1(3*m-2:3*m,3*t1(k-1,1)-2:3*t1(k-1,1))=sigma1(3*m-2:3*m,1:3);                % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))')
                        P1(:,3*t1(k-1,1)-2:3*t1(k-1,1))=P1(:,1:3); % P1=x(k)*x(k);
                        t1(k-1,1)=t1(k-1,1)-1;
                    else
                        tru1(3*m-2:3*m,t1(k-1,1)+i)=(A+F1*F(t1(k-1,1)+i-1,1)*E)*x2(3*m-2:3*m,t1(k-1,1)+i-1)+...
                            B*w1(t1(k-1,1)+i-1,1);
                        z1(t1(k-1,1)+i,m)=(C(m,:)+H1*F(t1(k-1,1)+i,1)*E)*tru1(3*m-2:3*m,t1(k-1,1)+i)+v1(t1(k-1,1)+i,m);
                        x2(3*m-2:3*m,t1(k-1,1)+i)=x1(3*m-2:3*m,t1(k-1,1)+i)+K1(3*m-2:3*m,t1(k-1,1)+i)*...
                            (z1(t1(k-1,1)+i,m)-C1(m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*x1(3*m-2:3*m,t1(k-1,1)+i));         % x2=x(k|k);
                        x1(3*m-2:3*m,t1(k-1,1)+i+1)=A1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*x1(3*m-2:3*m,t1(k-1,1)+i)+...
                            L1(3*m-2:3*m,t1(k-1,1)+i)*(z1(t1(k-1,1)+i,m)-...
                            C1(m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*x1(3*m-2:3*m,t1(k-1,1)+i));           % x1=x(k|k-1);
                    end
                else
                    x2(3*m-2:3*m,t1(k,1))=x1(3*m-2:3*m,t1(k,1))+K1(3*m-2:3*m,t1(k,1))*...
                        (z1(t1(k,1),m)-C1(m,3*t1(k,1)-2:3*t1(k,1))*x1(3*m-2:3*m,t1(k,1)));         % x2=x(k|k);
                    x1(3*m-2:3*m,t1(k,1)+1)=A1(3*m-2:3*m,3*t1(k,1)-2:3*t1(k,1))*x1(3*m-2:3*m,t1(k,1))+...
                        L1(3*m-2:3*m,t1(k,1))*(z1(t1(k,1),m)-...
                        C1(m,3*(t1(k,1))-2:3*(t1(k,1)))*x1(3*m-2:3*m,t1(k,1)));           % x1=x(k|k-1);
                end
                % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))')
                sigma1(3*m-2:3*m,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1))=...
                    A*sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*...
                    (eye(3)+E'*inv(M2(t1(k-1,1)+i,m))*E*sigma1(3*m-2:3*m,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*A'-...
                    delta1(3*m-2:3*m,t1(k-1,1)+i)*inv(Xi1(t1(k-1,1)+i,m))*delta1(3*m-2:3*m,t1(k-1,1)+i)'+...
                    B*Q_k*B'+(1/a)*F1*F1';  
                % P1=x(k)*x(k)'; 
                P1(:,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1))=...
                    A*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*A'+...
                    A*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(M1(t1(k-1,1)+i,1))*E*...
                    P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*A'+(1/a)*F1*F1'+B*Q_k*B';
            end
            if tau1(k,1)==0
                x3(3*m-2:3*m,k)=x2(3*m-2:3*m,t1(k,1));
            elseif tau1(k,1)==1
                x3(3*m-2:3*m,k)=x1(3*m-2:3*m,t1(k,1)+1);
            else
                x3(3*m-2:3*m,k)=((N-tau1(k,1)+1)/N)*x1(3*m-2:3*m,t1(k,1)+1);
            end
        end 
    end
    st(m,2)=cputime-st(m,1);
 %    TrP1(i,1)=trace(P2(3*i-2:3*i,1:3));
end

sigma2(1:9,1:9)=sigma1(1:9,1:9);
st(4,1)=cputime;
for k=1:iter-1
    for i=1:L
        for j=1:L
                sigma2(3*i-2:3*i,3*(j+L*k)-2:3*(j+L*k))=((N-tau1(k,1)+1)/N)^2.*...
                    (A-L1(3*i-2:3*i,k)*C(i,:))*sigma2(3*i-2:3*i,3*(j+L*(k-1))-2:3*(j+L*(k-1)))*(A-L1(3*j-2:3*j,k)*C(j,:))'+...
                    B*Q_k*B'+L1(3*i-2:3*i,k)*R_k(i,j)*L1(3*j-2:3*j,k)'+...
                    (A-A1(3*i-2:3*i,3*k-2:3*k)+L1(3*i-2:3*i,k)*(C1(i,3*k-2:3*k)-C(i,:)))*...
                    (P1(:,3*k-2:3*k)-sigma2(3*i-2:3*i,3*(j+L*(k-1))-2:3*(j+L*(k-1))))*...
                    (A-A1(3*j-2:3*j,3*k-2:3*k)+L1(3*j-2:3*j,k)*(C1(j,3*k-2:3*k)-C(j,:)))+...
                    ((A-L1(3*i-2:3*i,k)*C(i,:))*sigma2(3*i-2:3*i,3*(i+L*(k-1))-2:3*(i+L*(k-1)))+...
                    (A-A1(3*i-2:3*i,3*k-2:3*k)+L1(3*i-2:3*i,k)*(C1(i,3*k-2:3*k)-C(i,:)))*...
                    (P1(:,3*k-2:3*k)-sigma2(3*i-2:3*i,3*(i+L*(k-1))-2:3*(i+L*(k-1)))))*...
                    E'*inv(M2(k,i))*E*...                    
                    ((A-L1(3*i-2:3*i,k)*C(i,:))*sigma2(3*i-2:3*i,3*(i+L*(k-1))-2:3*(i+L*(k-1)))+...
                    (A-A1(3*i-2:3*i,3*k-2:3*k)+L1(3*i-2:3*i,k)*(C1(i,3*k-2:3*k)-C(i,:)))*...
                    (P1(:,3*k-2:3*k)-sigma2(3*i-2:3*i,3*(i+L*(k-1))-2:3*(i+L*(k-1)))))'+...
                    +(1/a)*(F1-L1(3*i-2:3*i,k)*H1)*(F1-L1(3*j-2:3*j,k)*H1)'-...
                    B*S_k(j,1)*L1(3*j-2:3*j,k)'-L1(3*i-2:3*i,k)*S_k(i,1)'*B';
        end
    end
      Pi1(:,3*k-2:3*k)=inv(I1'*inv(sigma2(:,9*k-8:9*k))*I1);
      if  isnan(Pi1(:,3*k-2:3*k))
          Pi1(:,3*k-2:3*k)=P;
      end
      Trtheta2(k,1)=trace(Pi1(:,3*k-2:3*k));             % Trace(P1);
      Omega1(:,9*k-8:9*k)=Pi1(:,3*k-2:3*k)*I1'*inv(sigma2(:,9*k-8:9*k));
      x4(:,k)=Omega1(:,9*k-8:9*k)*x3(:,k) ;
end
st(4,2)=cputime-st(4,1);
iter=300;
for i=1:L
    for k=1:iter
        fv(3*i-2:3*i,3*k-2:3*k)=0.1*(tru(:,k)-x3(3*i-2:3*i,k))*(tru(:,k)-x3(3*i-2:3*i,k))';
        b1(3*i-2,k)=fv(3*i-2,3*k-2);
        b1(3*i-1,k)=fv(3*i-1,3*k-1);
        b1(3*i,k)=fv(3*i,3*k);
        d1(3*i-2,1)=d1(3*i-2,1)+b1(3*i-2,k);
        d1(3*i-1,1)=d1(3*i-1,1)+b1(3*i-1,k);
        d1(3*i,1)=d1(3*i,1)+b1(3*i,k);
    end                                                                                                                                                                                                                                                                                                                                                                    
end
for k=1:iter
    b1(10,k)=Pi1(1,3*k-2);
    b1(11,k)=Pi1(2,3*k-1);
    b1(12,k)=Pi1(3,3*k);
    d1(10,1)=d1(10,1)+b1(10,k);
    d1(11,1)=d1(11,1)+b1(11,k);
    d1(12,1)=d1(12,1)+b1(12,k);        
end
for i=1:12
    d1(i,1)=d1(i,1)/iter;
end

%%%%%%%%%%%%%%% Position states %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(tru(1,1:iter),'r');       hold on;    plot(x3(1,1:iter),'b');       hold on; 
plot(x3(4,1:iter),'c');        hold on;    plot(x3(7,1:iter),'Color',[1 0.5 0]);       hold on; 
xlabel('t/step');     ylabel('States');
legend('Actual state of position','Estimated state of sensor 1',...
    'Estimated state of sensor 2','Estimated state of sensor 3');

%%%%%%%%%%%%%%% Velocity states %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(tru(2,1:iter),'r');       hold on;    plot(x3(2,1:iter),'b');       hold on; 
plot(x3(5,1:iter),'c');        hold on;    plot(x3(8,1:iter),'Color',[1 0.5 0]);       hold on;
xlabel('t/step');     ylabel('States');
legend('Actual state of velocity','Estimated state of sensor 1',...
    'Estimated state of sensor 2','Estimated state of sensor 3');
%%%%%%%%%%%%%%% Acceleration states %%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(tru(3,1:iter),'r');       hold on;    plot(x3(3,1:iter),'b');       hold on;
plot(x3(6,1:iter),'c');        hold on;    plot(x3(9,1:iter),'Color',[1 0.5 0]);       hold on;
xlabel('t/step');     ylabel('States');
legend('Actual state of acceleration','Estimated state of sensor 1',...
    'Estimated state of sensor 2','Estimated state of sensor 3');
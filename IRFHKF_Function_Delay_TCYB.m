% Robust Finite Horizon Kalman Filtering_2015_Information Science
function [sigma1, Trtheta1, P1, x2, MIS1, MIS2, fv, et]=IRFHKF_Function_Delay_TCYB(T, A, B, C, E, F1, H1, F, a, beta,...
    Q_k, R_k, S_k, M1, M2, tru, z, x1, sigma1, P1, iter, tau, N) 
st=cputime;        % start runtime

for k=1:iter 
  %  tau(k,1)=1;
    t(k,1)=k-tau(k,1);
end

for k=1:iter   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  ADD  %%%%%%%%%%%%%%%%
    if tau(k,1)==0 
        if t(k,1)==0
             t(k,1)=1;
        end
        MIS1(k,1)=M1(k,1);
        MIS2(k,1)=M2(k,1);
        lambda1(k,1)=C*sigma1(:,3*k-2:3*k)*(eye(3)+E'*inv(M2(k,1))*E*sigma1(:,3*k-2:3*k))*...
            C'+R_k;% lambda1=C*sigma1*(I+E'*M2*E*sigma1)*C'+R_k
        % Xi1=C*sigma1
        Xi1(k,1)=C*(eye(3)+E'*inv(M1(k,1))*E)*sigma1(:,3*k-2:3*k)*C'+R_k;
        % Filtering parameters K1, A1, L1 
        L1(:,k)=A*sigma1(:,3*k-2:3*k)*(eye(3)+E'*inv(M2(k,1))*E*sigma1(:,3*k-2:3*k))*...
            C'*inv(lambda1(k,1));
        A1(:,3*k-2:3*k)=A+(A-L1(:,k)*C)*sigma1(:,3*k-2:3*k)*E'*inv(M2(k,1))*E;    
        K1(:,k)=(eye(3)+sigma1(:,3*k-2:3*k)*E'*inv(M1(k,1))*E)*...
            sigma1(:,3*k-2:3*k)*C'*inv(Xi1(k,1));    
        x2(:,k)=x1(:,k)+K1(:,k)*(z(k,1)-C*x1(:,k));         % x2=x(k|k);
        x1(:,k+1)=A1(:,3*k-2:3*k)*x1(:,k)+L1(:,k)*(z(k,1)-C*x1(:,k));           % x1=x(k|k-1);
        Trtheta1(k,1)=trace(sigma1(:,3*k-2:3*k));
        % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))') 
        sigma1(:,3*(k+1)-2:3*(k+1))=A*sigma1(:,3*k-2:3*k)*...
            (eye(3)+E'*inv(M2(k,1))*E*sigma1(:,3*k-2:3*k))*A'-...
            sigma1(:,3*k-2:3*k)*(eye(3)+E'*inv(M2(k,1))*E*sigma1(:,3*k-2:3*k))*...
            C'*inv(lambda1(k,1))*C*sigma1(:,3*k-2:3*k)*...
            (eye(3)+E'*inv(M2(k,1))*E*sigma1(:,3*k-2:3*k))*A'+...
            B*Q_k*B'+(1/a)*F1*F1';
  %%%      b2(:,k+1)=eig(sigma1(:,3*(k+1)-2:3*(k+1)));
        % P1=x(k)*x(k)'; 
        P1(:,3*(k+1)-2:3*(k+1))=A*P1(:,3*k-2:3*k)*A'+...
            A*P1(:,3*k-2:3*k)*E'*inv(M1(k,1))*E*P1(:,3*k-2:3*k)*A'+...
            (1/a)*F1*F1'+B*Q_k*B';
 %%%%       b3(:,k+1)=eig(P1(:,3*(k+1)-2:3*(k+1)));
    elseif tau(k,1)==1
        if t(k,1)==0
            t(k,1)=1;
        end
        tru1(:,t(k,1))=tru(:,t(k,1));
        z1(t(k,1),1)=z(t(k,1),1);   % z1=(C+HFE)tru(k)+v_k
        MIS1(t(k,1),1)=(1/a)-E*P1(:,3*t(k,1)-2:3*t(k,1))*E';       % M1=1/a-E*P*E'
        MIS2(t(k,1),1)=(1/a)-E*sigma1(:,3*t(k,1)-2:3*t(k,1))*E';   % M2=1/a-E*sigma1*E'
        lambda1(t(k,1),1)=C*sigma1(:,3*t(k,1)-2:3*t(k,1))*(eye(3)+...
             E'*inv(MIS2(t(k,1),1))*E*sigma1(:,3*t(k,1)-2:3*t(k,1)))*C'+R_k;% lambda1=C*sigma1*(I+E'*M2*E*sigma1)*C'+R_k
      
        Xi1(t(k,1),1)=C*(eye(3)+E'*inv(MIS1(t(k,1),1))*E)*sigma1(:,3*t(k,1)-2:3*t(k,1))*C'+R_k;            
            
        L1(:,t(k,1))=A*sigma1(:,3*t(k,1)-2:3*t(k,1))*(eye(3)+E'*inv(MIS2(t(k,1),1))*E*sigma1(:,3*t(k,1)-2:3*t(k,1)))*...
            C'*inv(lambda1(t(k,1),1));
        A1(:,3*t(k,1)-2:3*t(k,1))=A+(A-L1(:,t(k,1))*C)*sigma1(:,3*t(k,1)-2:3*t(k,1))*E'*inv(MIS2(t(k,1),1))*E;    
        K1(:,t(k,1))=(eye(3)+sigma1(:,3*t(k,1)-2:3*t(k,1))*E'*inv(MIS1(t(k,1),1))*E)*...
            sigma1(:,3*t(k,1)-2:3*t(k,1))*C'*inv(Xi1(t(k,1),1));    
        x2(:,t(k,1))=x1(:,t(k,1))+K1(:,t(k,1))*(z(t(k,1),1)-C*x1(:,t(k,1)));         % x2=x(k|k);
        x1(:,t(k,1)+1)=A1(:,3*t(k,1)-2:3*t(k,1))*x1(:,t(k,1))+L1(:,t(k,1))*(z(t(k,1),1)-C*x1(:,t(k,1)));           % x1=x(k|k-1);
        Trtheta1(t(k,1),1)=trace(sigma1(:,3*t(k,1)-2:3*t(k,1)));
        sigma1(:,3*(t(k,1)+1)-2:3*(t(k,1)+1))=A*sigma1(:,3*t(k,1)-2:3*t(k,1))*...
            (eye(3)+E'*inv(MIS2(t(k,1),1))*E*sigma1(:,3*t(k,1)-2:3*t(k,1)))*A'-...
            sigma1(:,3*t(k,1)-2:3*t(k,1))*(eye(3)+E'*inv(MIS2(t(k,1),1))*E*sigma1(:,3*t(k,1)-2:3*t(k,1)))*...
            C'*inv(lambda1(t(k,1),1))*C*sigma1(:,3*t(k,1)-2:3*t(k,1))*...
            (eye(3)+E'*inv(MIS2(t(k,1),1))*E*sigma1(:,3*t(k,1)-2:3*t(k,1)))*A'+...
            B*Q_k*B'+(1/a)*F1*F1';
         P1(:,3*(t(k,1)+1)-2:3*(t(k,1)+1))=A*P1(:,3*t(k,1)-2:3*t(k,1))*A'+...
            A*P1(:,3*t(k,1)-2:3*t(k,1))*E'*inv(MIS1(t(k,1),1))*E*P1(:,3*t(k,1)-2:3*t(k,1))*A'+...
            (1/a)*F1*F1'+B*Q_k*B';
    elseif tau(k,1)>1
        if t(k,1)==0
            t(k,1)=1;
        end
        t1(k,1)=t(k,1);
        tru1(:,t1(k,1))=tru(:,t1(k,1));
        z1(t1(k,1),1)=z(t1(k,1),1);   % z1=(C+HFE)tru(k)+v_k
        for i=1:tau(k,1)
            MIS1(t1(k-1,1)+i+1,1)=(1/a)-E*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E';       % M1=1/a-E*P*E'
            MIS2(t1(k-1,1)+i+1,1)=(1/a)-E*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E';   % M2=1/a-E*sigma1*E'
            lambda1(t1(k-1,1)+i+1,1)=C*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*(eye(3)+...
                E'*inv(MIS2(t1(k-1,1)+i,1))*E*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*C'+R_k;% lambda1=C*sigma1*(I+E'*M2*E*sigma1)*C'+R_k   
% Original:  Xi1(k,1)=C*(eye(3)+E'*inv(M1(k,1))*E)*sigma1(:,3*k-2:3*k)*C'+R_k;
            Xi1(t1(k-1,1)+i+1,1)=C*(eye(3)+E'*inv(MIS1(t1(k-1,1)+i,1))*E)*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*C'+R_k;
 % Filtering parameters K1, A1, L1 
            L1(:,t1(k-1,1)+i+1)=A*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*(eye(3)+...
                E'*inv(MIS2(t1(k-1,1)+i,1))*E*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*C'*inv(lambda1((t1(k-1,1)+i),1));
            A1(:,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1))=A+(A-L1(:,t1(k-1,1)+i)*C)*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(MIS2((t1(k-1,1)+i),1))*E;
            K1(:,t1(k-1,1)+i+1)=(eye(3)+sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(MIS1((t1(k-1,1)+i),1))*E)*...
                sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*C'*inv(Xi1((t1(k-1,1)+i),1)); 
            x2(:,t1(k-1,1)+i+1)=x1(:,t1(k-1,1)+i)+K1(:,(t1(k-1,1)+i))*(z(t1(k-1,1)+i,1)-C*x1(:,t1(k-1,1)+i));         % x2=x(k|k);
            x1(:,t1(k-1,1)+i+1)=A1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*x1(:,t1(k-1,1)+i)+L1(:,t1(k-1,1)+i)*(z(t1(k-1,1)+i,1)-C*x1(:,t1(k-1,1)+i));           % x1=x(k|k-1);
            Trtheta1(t1(k-1,1)+i+1,1)=trace(sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)));
        % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))') 
            sigma1(:,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1))=A*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*...
                (eye(3)+E'*inv(MIS2(t1(k-1,1)+i,1))*E*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*A'-...
                sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*(eye(3)+E'*inv(MIS2(t1(k-1,1)+i,1))*E*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*...
                C'*inv(lambda1(t1(k-1,1)+i,1))*C*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*...
                (eye(3)+E'*inv(MIS2(t1(k-1,1)+i,1))*E*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*A'+...
                B*Q_k*B'+(1/a)*F1*F1';  
            % P1=x(k)*x(k)'; 
            P1(:,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1))=A*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*A'+...
                A*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(MIS1(t1(k-1,1)+i,1))*E*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*A'+...
                (1/a)*F1*F1'+B*Q_k*B';
        end
    end
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  End ADD  %%%%%%%%%%%%
for k=1:iter
    fv(:,3*k-2:3*k)=0.1*(tru(:,k)-x2(:,k))*(tru(:,k)-x2(:,k))';
end
et=cputime-st;        % runtime

% A distributed fusion estimation approach based on the robust finite horizon filtering without packet disorders
function [sigma1, Trtheta1, P1, x3, M1, M2, fv, et]=RFHDFE_LZOH_Function_TCYB(T, A, B, C, E, F1, H1, F, a, beta,...
    w1, v1, Q_k, R_k, S_k, tru, z, x1, x2, sigma1, P1, iter, tau, N) 
st=cputime;        % start runtime
for k=1:iter
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
end

for k=1:iter
    if t1(k,1)==0
        x2(:,t1(k,1)+1)=tru(:,1);
        x3(:,k)=tru(:,1);
        if tau1(k,1)>=1
            x3(:,k)=((N-tau1(k,1)+1)/N)*tru(:,1);
        end
        continue;
    else
        tru1(:,t1(k,1))=tru(:,t1(k,1));
        z1(t1(k,1),1)=z(t1(k,1),1);   % z1=(C+HFE)tru(k)+v_k
    end      
    if s1(k,1)<=1
        if s1(k,1)==1
            M1(t1(k,1),1)=(1/a)-E*P1(:,3*t1(k,1)-2:3*t1(k,1))*E';       % M1=1/a-E*P*E'
            M2(t1(k,1),1)=(1/a)-E*sigma1(:,3*t1(k,1)-2:3*t1(k,1))*E';   % M2=1/a-E*sigma1*E'
            lambda1(:,t1(k,1))=(eye(3)+sigma1(:,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M2(t1(k,1),1))*E)*...
                sigma1(:,3*t1(k,1)-2:3*t1(k,1))*C';% lambda1=sigma1*C'
            delta1(:,t1(k,1))=A*sigma1(:,3*t1(k,1)-2:3*t1(k,1))*...
                (eye(3)+E'*inv(M2(t1(k,1),1))*E*sigma1(:,3*t1(k,1)-2:3*t1(k,1)))*C'+...
                (1/a)*F1*H1'+B*S_k;% delta1=A*sigma1
            % Xi1=C*sigma1
            Xi1(t1(k,1),1)=C*sigma1(:,3*t1(k,1)-2:3*t1(k,1))*...
                (eye(3)+E'*inv(M2(t1(k,1),1))*E*sigma1(:,3*t1(k,1)-2:3*t1(k,1)))*C'+...
                (1/a)*H1*H1'+R_k;
            % theta1=E((x(k)-x(k|k))(x(k)-x(k|k))')
            theta1(:,3*t1(k,1)-2:3*t1(k,1))=sigma1(:,3*t1(k,1)-2:3*t1(k,1))+...
                sigma1(:,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M1(t1(k,1),1))*E*sigma1(:,3*t1(k,1)-2:3*t1(k,1))-...
                lambda1(:,t1(k,1))*inv(Xi1(t1(k,1),1))*lambda1(:,t1(k,1))';
            b1(:,t1(k,1))=eig(theta1(:,3*t1(k,1)-2:3*t1(k,1)));
            Trtheta1(t1(k,1),1)=trace(sigma1(:,3*t1(k,1)-2:3*t1(k,1)));             % Trace(P1);
            % Filtering parameters C1, K1, A1, L1 
            C1(1,3*t1(k,1)-2:3*t1(k,1))=C*(eye(3)+...
                sigma1(:,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M2(t1(k,1),1))*E);
            K1(:,t1(k,1))=lambda1(:,t1(k,1))*inv(Xi1(t1(k,1),1));    
            A1(:,3*t1(k,1)-2:3*t1(k,1))=A*(eye(3)+...
                sigma1(:,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M2(t1(k,1),1))*E);
            L1(:,t1(k,1))=delta1(:,t1(k,1))*inv(Xi1(t1(k,1),1));
            x2(:,t1(k,1))=x1(:,t1(k,1))+K1(:,t1(k,1))*...
                (z1(t1(k,1),1)-C1(1,3*t1(k,1)-2:3*t1(k,1))*x1(:,t1(k,1)));         % x2=x(k|k);
            x1(:,t1(k,1)+1)=A1(:,3*t1(k,1)-2:3*t1(k,1))*x1(:,t1(k,1))+L1(:,t1(k,1))*...
                (z1(t1(k,1),1)-C1(1,3*t1(k,1)-2:3*t1(k,1))*x1(:,t1(k,1)));           % x1=x(k|k-1);
            % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))') 
            sigma1(:,3*(t1(k,1)+1)-2:3*(t1(k,1)+1))=A*sigma1(:,3*t1(k,1)-2:3*t1(k,1))*...
                (eye(3)+E'*inv(M2(t1(k,1),1))*E*sigma1(:,3*t1(k,1)-2:3*t1(k,1)))*A'-...
                delta1(:,t1(k,1))*inv(Xi1(t1(k,1),1))*delta1(:,t1(k,1))'+...
                B*Q_k*B'+(1/a)*F1*F1';  
            b2(:,t1(k,1)+1)=eig(sigma1(:,3*(t1(k,1)+1)-2:3*(t1(k,1)+1)));
            % P1=x(k)*x(k)'; 
            P1(:,3*(t1(k,1)+1)-2:3*(t1(k,1)+1))=A*P1(:,3*t1(k,1)-2:3*t1(k,1))*A'+...
                A*P1(:,3*t1(k,1)-2:3*t1(k,1))*E'*inv(M1(t1(k,1),1))*E*P1(:,3*t1(k,1)-2:3*t1(k,1))*A'+...
                (1/a)*F1*F1'+B*Q_k*B';
            b3(:,t1(k,1)+1)=eig(P1(:,3*(t1(k,1)+1)-2:3*(t1(k,1)+1)));
        end
        if tau1(k,1)==0
            x3(:,k)=x2(:,t1(k,1));
        elseif tau1(k,1)==1
            x3(:,k)=x1(:,t1(k,1)+1);
        else
            x3(:,k)=((N-tau1(k,1)+1)/N)*x1(:,t1(k,1)+1);
        end
    end
    % if s1(k,1)>1 || tau1(k,1)>1
    if s1(k,1)>1
        for i=1:s1(k,1)
            M1(t1(k-1,1)+i,1)=(1/a)-E*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E';       % M1=1/a-E*P*E'
            M2(t1(k-1,1)+i,1)=(1/a)-E*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E';   % M2=1/a-E*sigma1*E'
            lambda1(:,t1(k-1,1)+i)=(eye(3)+sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*...
                inv(M2(t1(k-1,1)+i,1))*E)*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*C';% lambda1=sigma1*C'
            delta1(:,t1(k-1,1)+i)=A*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*...
                (eye(3)+E'*inv(M2(t1(k-1,1)+i,1))*E*...
                sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*C'+(1/a)*F1*H1'+B*S_k;% delta1=A*sigma1
            % Xi1=C*sigma1
            Xi1(t1(k-1,1)+i,1)=C*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*...
                (eye(3)+E'*inv(M2(t1(k-1,1)+i,1))*E*...
                sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*C'+(1/a)*H1*H1'+R_k;
            % theta1=E((x(k)-x(k|k))(x(k)-x(k|k))')
            theta1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))=sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))+...
                sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(M1(t1(k-1,1)+i,1))*E*...
                sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))-...
                lambda1(:,t1(k-1,1)+i)*inv(Xi1(t1(k-1,1)+i,1))*lambda1(:,t1(k-1,1)+i)';
            b1(:,t1(k-1,1)+i)=eig(theta1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)));
            Trtheta1(t1(k-1,1)+i,1)=trace(theta1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)));% Trace(P1);
            % Filtering parameters C1, K1, A1, L1 
            C1(1,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))=C*(eye(3)+...
                sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(M2(t1(k-1,1)+i,1))*E);
            K1(:,t1(k-1,1)+i)=lambda1(:,t1(k-1,1)+i)*inv(Xi1(t1(k-1,1)+i,1));    
            A1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))=A*(eye(3)+...
                sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(M2(t1(k-1,1)+i,1))*E);
            L1(:,t1(k-1,1)+i)=delta1(:,t1(k-1,1)+i)*inv(Xi1(t1(k-1,1)+i,1));
            if i<s1(k,1)
                if t1(k-1,1)==0
                    t1(k-1,1)=t1(k-1,1)+1;
                    tru1(:,t1(k-1,1))=tru(:,1);
                    z1(t1(k-1,1),1)=z(1,1);
                    x2(:,t1(k-1,1))=x1(:,t1(k-1,1))+K1(:,t1(k-1,1))*...
                        (z1(t1(k-1,1),1)-C1(1,3*t1(k-1,1)-2:3*t1(k-1,1))*x1(:,t1(k-1,1)));         % x2=x(k|k);
                    x1(:,t1(k-1,1)+1)=A1(:,3*t1(k-1,1)-2:3*t1(k-1,1))*x1(:,t1(k-1,1))+...
                        L1(:,t1(k-1,1))*(z1(t1(k-1,1),1)-...
                        C1(1,3*t1(k-1,1)-2:3*t1(k-1,1))*x1(:,t1(k-1,1)));           % x1=x(k|k-1);
                    sigma1(:,3*t1(k-1,1)-2:3*t1(k-1,1))=sigma1(:,1:3);                % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))')
                    P1(:,3*t1(k-1,1)-2:3*t1(k-1,1))=P1(:,1:3); % P1=x(k)*x(k);
                    t1(k-1,1)=t1(k-1,1)-1;
                else
                    tru1(:,t1(k-1,1)+i)=(A+F1*F(t1(k-1,1)+i-1,1)*E)*x2(:,t1(k-1,1)+i-1)+...
                        B*w1(t1(k-1,1)+i-1,1);
                    z1(t1(k-1,1)+i,1)=(C+H1*F(t1(k-1,1)+i,1)*E)*tru1(:,t1(k-1,1)+i)+v1(t1(k-1,1)+i,1);
                    x2(:,t1(k-1,1)+i)=x1(:,t1(k-1,1)+i)+K1(:,t1(k-1,1)+i)*...
                        (z1(t1(k-1,1)+i,1)-C1(1,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*x1(:,t1(k-1,1)+i));         % x2=x(k|k);
                    x1(:,t1(k-1,1)+i+1)=A1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*x1(:,t1(k-1,1)+i)+...
                        L1(:,t1(k-1,1)+i)*(z1(t1(k-1,1)+i,1)-...
                        C1(1,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*x1(:,t1(k-1,1)+i));           % x1=x(k|k-1);
                end
            else
                x2(:,t1(k,1))=x1(:,t1(k,1))+K1(:,t1(k,1))*...
                    (z1(t1(k,1),1)-C1(1,3*t1(k,1)-2:3*t1(k,1))*x1(:,t1(k,1)));         % x2=x(k|k);
                x1(:,t1(k,1)+1)=A1(:,3*t1(k,1)-2:3*t1(k,1))*x1(:,t1(k,1))+...
                    L1(:,t1(k,1))*(z1(t1(k,1),1)-...
                    C1(1,3*(t1(k,1))-2:3*(t1(k,1)))*x1(:,t1(k,1)));           % x1=x(k|k-1);
            end
            % sigma1=E((x(k)-x(k|k-1))(x(k)-x(k|k-1))')
            sigma1(:,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1))=...
                A*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*...
                (eye(3)+E'*inv(M2(t1(k-1,1)+i,1))*E*sigma1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i)))*A'-...
                delta1(:,t1(k-1,1)+i)*inv(Xi1(t1(k-1,1)+i,1))*delta1(:,t1(k-1,1)+i)'+...
                B*Q_k*B'+(1/a)*F1*F1';  
            b2(:,t1(k-1,1)+i+1)=eig(sigma1(:,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1)));
            % P1=x(k)*x(k)'; 
            P1(:,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1))=...
                A*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*A'+...
                A*P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*E'*inv(M1(t1(k-1,1)+i,1))*E*...
                P1(:,3*(t1(k-1,1)+i)-2:3*(t1(k-1,1)+i))*A'+(1/a)*F1*F1'+B*Q_k*B';
            b3(:,t1(k-1,1)+i+1)=eig(P1(:,3*(t1(k-1,1)+i+1)-2:3*(t1(k-1,1)+i+1)));
        end
        if tau1(k,1)==0
            x3(:,k)=x2(:,t1(k,1));
        elseif tau1(k,1)==1
            x3(:,k)=x1(:,t1(k,1)+1);
        else
            x3(:,k)=((N-tau1(k,1)+1)/N)*x1(:,t1(k,1)+1);
        end
    end 
end
for k=1:iter
    fv(:,3*k-2:3*k)=0.1*(tru(:,k)-x3(:,k))*(tru(:,k)-x3(:,k))';
end
et=cputime-st;        % runtime
end
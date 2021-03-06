% Robust Finite Horizon Kalman Filtering
function [sigma1, Trtheta1, P1, x2]=IRFHKF_Function(T, A, B, C, E, F1, H1, F, a, beta,...
    Q_k, R_k, S_k, M1, M2, tru, z, x1, sigma1, P1, iter) 
for k=1:iter   
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
    b2(:,k+1)=eig(sigma1(:,3*(k+1)-2:3*(k+1)));
    % P1=x(k)*x(k)'; 
    P1(:,3*(k+1)-2:3*(k+1))=A*P1(:,3*k-2:3*k)*A'+...
        A*P1(:,3*k-2:3*k)*E'*inv(M1(k,1))*E*P1(:,3*k-2:3*k)*A'+...
        (1/a)*F1*F1'+B*Q_k*B';
    b3(:,k+1)=eig(P1(:,3*(k+1)-2:3*(k+1)));
end

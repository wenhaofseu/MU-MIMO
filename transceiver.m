%Statistically Robust Transceiver Design for Multi-RIS Assisted Multi-User MIMO Systems
function [G,F,sumMSE_all] = transceiver(H_eff_est,F,Pt,sigma_q,iter_max,eplision)
    [Nk,M,K] = size(H_eff_est);
    [~,Dk] = size(F);
    Dk = Dk/K;
    sumMSE_all = zeros(iter_max,1);
    for iter = 1:iter_max
        %% update G and calculate C
        G = zeros(Dk,Nk,K);
        C = zeros(M);
        for k = 1:K
            Fk = F(:,1+(k-1)*Dk:k*Dk);
            Hk_eff_est = H_eff_est(:,:,k);
            Ak = (Hk_eff_est*F)*(Hk_eff_est*F)'+sigma_q*eye(Nk);
            Gk = Fk'*Hk_eff_est'/Ak;
            G(:,:,k) = Gk;
            C = C + (Gk*Hk_eff_est)'*(Gk*Hk_eff_est);
        end
        [Uc,LAMBDAc]=eig(C);
        B = Uc'*C*Uc;
        %% check lambda = 0 is optimal?
        flag = 0;
        if rank(C)==M
            temp1 = 0;
            for m = 1:M
                temp1 = temp1+B(m,m)/LAMBDAc(m,m)^2;
            end
            if temp1<=Pt
                flag = 1;
            end
        end
        %% calculate lambda
        lambda_min = 0;
        lambda_max = sqrt(trace(B)/Pt);
        while(flag==0)
            lambda = (lambda_min+lambda_max)/2;
            temp2 = 0;
            for m = 1:M
                temp2 = temp2+B(m,m)/(LAMBDAc(m,m)+lambda)^2;
            end
            temp2 = real(temp2);
            if(temp2<Pt)
                lambda_max = lambda;
            else
                lambda_min = lambda;
            end
            %% bisection converage
            if abs(lambda_max-lambda_min)<1e-3 && abs(temp2-Pt)<1e-5
                break
            elseif abs(lambda_max-lambda_min)<1e-7 && lambda_min==0
                break
            end
        end
        %% update F
        F = zeros(M,Dk*K);
        for k = 1:K
            Gk = G(:,:,k);
            Hk_eff_est = H_eff_est(:,:,k);
            Fk = (C+lambda*eye(M))\(Gk*Hk_eff_est)';
            F(:,1+(k-1)*Dk:k*Dk) = Fk;
        end
        norm(F,'fro')^2
        %% calculate MSE
        MSE = zeros(K,1);
        for k = 1:K
            Gk = G(:,:,k);
            Fk = F(:,1+(k-1)*Dk:k*Dk);
            Hk_eff_est = H_eff_est(:,:,k);
            MSEk = trace((Gk*Hk_eff_est*F)*(Gk*Hk_eff_est*F)')-2*real(trace(Gk*Hk_eff_est*Fk))+Dk+sigma_q*trace(Gk*Gk');
            MSE(k) = MSEk;
        end
        sumMSE = sum(MSE);
        sumMSE_all(iter) = sumMSE;
        %% algorithm converage
        if iter>=5 && abs(sumMSE-sumMSE_all(iter-1))/abs(sumMSE)<eplision
            sumMSE_all = sumMSE_all(1:iter);
            break
        end
    end
end


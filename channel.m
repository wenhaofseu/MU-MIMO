%Statistically Robust Transceiver Design for Multi-RIS Assisted
%Multi-user MIMO Systems
clc;clear;close all;
rng(1);
M = 32;%number of BS antennas
Lr = 64;%number of RISr elements
R = 2;
K = 4;%number of users
Nk = 4;%number of userk antennas
BS = [0,0,10];%location of BS(m)
RIS1 = [-10,50,10];%location of RIS1(m)
RIS2 = [10,50,10];%location of RIS2(m)
radius = 5;
center = [0,50,1.5];
Users = sample(radius,center,K);%location of Users(m)
rho = 0.5;
%% model:H = Phi^(1/2)*Hw*Psai^(1/2)
%% normalized correlation matrix
Psai_dk = correlation(rho,M);
Phi_dk = correlation(rho,Nk);
Psai_b = correlation(rho,M);
Phi_br = correlation(rho,Lr);
Psai_rk = correlation(rho,Lr);
Phi_rk = correlation(rho,Nk);
%% use blkdiag
Phi_b = zeros(Lr*R);
for r = 1:R
    Phi_b(1+(r-1)*Lr:r*Lr,1+(r-1)*Lr:r*Lr) = correlation(rho,Lr);
end
Psai_k = zeros(Lr*R);
for r = 1:R
    Psai_k(1+(r-1)*Lr:r*Lr,1+(r-1)*Lr:r*Lr) = correlation(rho,Lr);
end
%% path loss
beta0 = 30;%dB
d0 = 1;%(m)
gamma_dk = 3.6;
gamma_br = 2.5;
gamma_rk = 2;
d = @(x1,x2) norm(x1-x2);%calculate distance
beta = @(d,gamma) 10^(-beta0/10)*(d/d0)^(-gamma);%calculate path loss(linear)
%% generate beta_dk
beta_d = zeros(K,1);
for k = 1:K
    userk = Users(k,:);
    beta_dk = beta(d(BS,userk),gamma_dk);
    beta_d(k) = beta_dk;
end
%% generate beta_br
beta_b = zeros(1,R);
for r = 1:R
    if r == 1
        RIS = RIS1;
    else
        RIS = RIS2;
    end
    beta_br = beta(d(BS,RIS),gamma_br);
    beta_b(r) = beta_br;
end
%% generate beta_rk
beta_RK = zeros(R,K);
for r = 1:R
    if r == 1
        RIS = RIS1;
    else
        RIS = RIS2;
    end
    for k = 1:K
        userk = Users(k,:);
        beta_rk = beta(d(RIS,userk),gamma_rk);
        beta_RK(r,k) = beta_rk;
    end
end
number = 100;%number of channel realization
%% store channel matrix
Hd_est_all = zeros(Nk,M,K,number);
Hb_est_all = zeros(R*Lr,M,number);
HRK_est_all = zeros(Nk,R*Lr,K,number);
for num = 1:number
    for k = 1:K
        %% generate Hdk_est
        beta_dk = beta_d(k);
        Hdk_w = sqrt(1/2)*(randn(Nk,M)+1i*randn(Nk,M));
        Hdk_est = sqrt(beta_dk)*Phi_dk^(1/2)*Hdk_w*Psai_dk^(1/2);
        Hd_est_all(:,:,k,num) = Hdk_est;
    end
    for r = 1:R
        %% generate Hbr_est
        beta_br = beta_b(r);
        Hbr_w = sqrt(1/2)*(randn(Lr,M)+1i*randn(Lr,M));
        Hbr_est = sqrt(beta_br)*Phi_br^(1/2)*Hbr_w*Psai_b^(1/2);
        Hb_est_all(1+(r-1)*Lr:r*Lr,:,num) = Hbr_est;
    end
    for r = 1:R
        for k = 1:K
            %% generate Hrk_est
            beta_rk = beta_RK(r,k);
            Hrk_w = sqrt(1/2)*(randn(Nk,Lr)+1i*randn(Nk,Lr));
            Hrk_est = sqrt(beta_rk)*Phi_rk^(1/2)*Hrk_w*Psai_rk^(1/2);
            HRK_est_all(:,1+(r-1)*Lr:r*Lr,k,num) = Hrk_est;
        end
    end   
end
save('channel.mat','Hd_est_all','Hb_est_all','HRK_est_all');
clc;clear;close all;
load('channel.mat','Hd_est_all','Hb_est_all','HRK_est_all');
M = 32;%number of BS
Lr = 64;%number of RISr elements
R = 2;%number of RIS
K = 4;%number of users
Nk = 4;%number of userk antennas
Dk = 2;%number of userk data streams
% sigma_q = 1e-8;%-110dbm
sigma_q = 1e-12;
Pt = 1;%0dB
Hd_est = Hd_est_all(:,:,:,1);
Hb_est = Hb_est_all(:,:,1);
HRK_est = HRK_est_all(:,:,:,1);
Theta = ones(Lr*R);%Phase-shifters
% %% set NMSE
% NMSE = 0;
%% define channel
H_eff_est = zeros(Nk,M,K);
for k = 1:K
    Hdk_est = Hd_est(:,:,k);
    Hk_est = HRK_est(:,:,k);
    Hk_eff_est = Hdk_est+Hk_est*Theta*Hb_est;
    H_eff_est(:,:,k) = Hk_eff_est;
end
iter_max = 200;
eplision = 1e-3;
for num = 1:100
    num
    %% initial F
    F = randn(M,Dk*K)+1i*randn(M,Dk*K);
    F = sqrt(Pt)*F/norm(F,'fro');
    [~,~,sumMSE_all] = transceiver(H_eff_est,F,Pt,sigma_q,iter_max,eplision);
    xlabel('iteration number');
    ylabel('Sum MSE of users');
    plot(1:length(sumMSE_all),sumMSE_all);hold on;
end
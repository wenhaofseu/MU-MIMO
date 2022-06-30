%"Zero-Forcing Methods for Downlink Spatial
%Multiplexing in Multiuser MIMO Channels"
clc;clear;close all;
Pt = 10;
sigma_n_q = 1;
nT = 16;
nRj = 4;
K = 4;
nR = K*nRj;
%% generate channel
Hs = [];
for j = 1:K
    Hj = randn(nRj,nT)+1i*randn(nRj,nT);
    Hs = [Hs;Hj];
end
%% Sum Capacity Block Diagonalization
Sigma = [];
D = 0;
Ms = [];
for j = 1:K
    Hj_tilde = [Hs(1:(j-1)*nRj,:);Hs(1+j*nRj:end,:)];
    Lj_tilde = rank(Hj_tilde);
    [Uj_tilde,Sigmaj_tilde,Vj_tilde] = svd(Hj_tilde);
    Vj_tilde0 = Vj_tilde(:,Lj_tilde+1:end);
    Hj = Hs(1+(j-1)*nRj:j*nRj,:);
    Lj_bar = rank(Hj*Vj_tilde0);
    D = D + Lj_bar;
    [Uj,Sigmaj,Vj] = svd(Hj*Vj_tilde0);
    Sigmaj = Sigmaj(1:Lj_bar,1:Lj_bar);
    Sigma = blkdiag(Sigma,Sigmaj);
    Vj1 = Vj(:,1:Lj_bar);
    Ms = [Ms,Vj_tilde0*Vj1];
end
[Sigma0,Index] = sort(diag(Sigma),'descend');
Sigma0 = diag(Sigma0);
L = D;
temp1 = sigma_n_q./diag(Sigma0^2);
while(1)
    level = (Pt+sum(temp1(1:L)))/L;%level = 1/p0
    if level>=temp1(L)
        break
    end
    L = L - 1;
end
p0 = [level - temp1(1:L);zeros(D-L,1)];
p = p0(Index);
Lambda = diag(p);
Ms = Ms*Lambda^(1/2);
C_BD1 = log2(det(eye(D)+Sigma^2*Lambda/sigma_n_q));
C_BD2 = log2(det(eye(nR)+(Hs*Ms)*(Hs*Ms)'/sigma_n_q));
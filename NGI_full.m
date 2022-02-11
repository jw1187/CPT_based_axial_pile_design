function [Q,Qf,Qp,Qt,Qz,Qzt,tau,qbz]=NGI_full(z,qc,Ic,sg,sgp,D,WT,L)

dz = z(2)-z(1);
Ar = 1-((D-2*WT)^2./D^2);
R = sqrt(D.^2-(D-2*WT)^2)/2;
Fsig = (sgp/100).^0.25;Fsig(isnan(Fsig)) = 0;
Dr = 0.4*log(qc./(22.*sqrt(sgp.*101.3)));Dr(isnan(Dr)) = 0;
Dr(find(Dr<0.1 & Dr~=0)) = 0.1;
y = ((ones(size(z))*L-z*ones(size(L)))./R);
cdg = 0.96;     % cyclic degradation factor
Fdr = 2.1*(Dr-0.1).^1.7;
S_ind = find(Ic<2.6);
C_ind = find(Ic>=2.6);
tau = zeros(size(y));
ten = zeros(size(y));
alpha = zeros(size(Ic));
Qz = zeros(size(y));
Qzt = zeros(size(y));
qbz = zeros(size(y));
S_prc = zeros(size(L'));
S_fc_prc = zeros(size(L'));
S_ft_prc = zeros(size(L'));
% sands
tau(S_ind,:) = 1.3.*(z(S_ind)*ones(size(L)))./(ones(size(z(S_ind)))*L).*101.3.*((Fsig(S_ind).*Fdr(S_ind))*ones(size(L)));
ten(S_ind,:) = (z(S_ind)*ones(size(L)))./(ones(size(z(S_ind)))*L).*101.3.*((Fsig(S_ind).*Fdr(S_ind))*ones(size(L)));
for i = 1:length(L)
    tau(S_ind,i) = 0.1*sgp(S_ind).*(tau(S_ind,i)<0.1*sgp(S_ind)) + tau(S_ind,i).*(tau(S_ind,i)>=0.1*sgp(S_ind));
    ten(S_ind,i) = 0.1*sgp(S_ind).*(ten(S_ind,i)<0.1*sgp(S_ind)) + ten(S_ind,i).*(ten(S_ind,i)>=0.1*sgp(S_ind));
end
% silts or clay
su = (qc-sg)/15; su(find(su<0)) = 0;%su±j­¢=0
sai = su./sgp;
ind1 = find(sai(C_ind)<=0.25);alpha(C_ind(ind1)) = 0.32*(14.4-10).^0.3;
ind2 = find(sai(C_ind)>1);alpha(C_ind(ind2)) = 0.5*(sai(C_ind(ind2))).^-0.3;
ind3 = find(sai(C_ind)<1 & sai(C_ind)>0.25);
alpha(C_ind(ind3)) = interp1([0.25,1],[0.32*(14.4-10)^0.3,0.5],sai(C_ind(ind3)));
alpha(isnan(su(C_ind)./sgp(C_ind)) | isinf(su(C_ind)./sgp(C_ind)))=0;sai(isnan(sai) | isinf(sai)) = 0 ;
tau(C_ind,:) = (alpha(C_ind).*sai(C_ind).*sgp(C_ind))*ones(size(L));
ten(C_ind,:) = tau(C_ind,:);
tau = tau.*cdg;
ten = ten.*cdg;
% bearing capacity
[~,ind] = min(abs(z*ones(size(L))-ones(size(z))*L));

for i = 1:length(ind)
    temp_tau(i,:) = sum(3*tau(1:ind(i),i).*dz);
    Qf(i,:) = sum(pi*D*dz*tau(1:ind(i),i));
    Qt(i,:) = sum(pi*D*dz*ten(1:ind(i),i));
    temp_qc = min([0.7*qc./(1+3*Dr.^2),qc*Ar+12*cumsum(tau(:,i))*dz./(D-2*WT).*(1-Ar)]')';
    temp_qc = temp_qc.*(Ic<2.6) + (Ic>=2.6).*9.*su;
    temp_qc(isnan(temp_qc)) = 0;
    Qz(1:ind(i),i) = cumsum(pi*D*dz*tau(1:ind(i),i)) + temp_qc(1:ind(i))*pi*(D.^2)/4;
    Qzt(1:ind(i),i) = cumsum(pi*D*dz*ten(1:ind(i),i));
    qbz(1:ind(i),i) =  temp_qc(1:ind(i));
    S_prc(i,:) = sum(S_ind<=ind(i))./ind(i);
    S_fc_prc(i,:) = sum(pi*D*dz*(S_ind<=ind(i)).*tau(S_ind,i))/Qf(i,:);
    S_ft_prc(i,:) = sum(pi*D*dz*(S_ind<=ind(i)).*ten(S_ind,i))/Qt(i,:);
end
Qp = min([0.7*qc(ind)./(1+3*Dr(ind).^2),qc(ind)*Ar+4*temp_tau./(D-2*WT).*(1-Ar)]')'*pi*D^2/4;
Qp = Qp.*(Ic(ind)<2.6) + (Ic(ind)>=2.6).*9.*su(ind)*pi*D^2/4;
Q = Qf+Qp;
end


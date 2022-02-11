function [Q,Qf,Qp,Qt,Qz,Qzt,tau,qbz,sgrc,sgrd,S_prc,S_fc_prc,S_ft_prc]=ICP05_full(z,qc,Ic,sgp,D,WT,L,st)

dz = z(2)-z(1);
K = round(3*D/dz);
R = sqrt(D.^2-(D-2*WT)^2)/2;
cdg = 0.96;     % cyclic degradation factor
y = ((ones(size(z))*L-z*ones(size(L)))./R);
y(find(y<=8)) = 8;
Ar = 1-((D-2*WT)^2./D^2);
S_ind = find(Ic<2.6);
C_ind = find(Ic>=2.6);
sgrc = zeros(size(y));
sgrd = zeros(size(y));
tau = zeros(size(y));
ten = zeros(size(y));
Qz = zeros(size(y));
Qzt = zeros(size(y));
qbz = zeros(size(y));
S_prc = zeros(size(L'));
S_fc_prc = zeros(size(L'));
S_ft_prc = zeros(size(L'));
% reduction_f = (Ic<2.05).*0.8 + (Ic>=2.05 & Ic<2.95).*0.67 + (Ic>=2.95).*1;  %fo ICP
delta_cv = (Ic<2.05).*30 + (Ic>=2.05 & Ic<2.6).*27 + (Ic>=2.6 & Ic<2.95).*25 + (Ic>=2.95).*20;
YSR = ones(size(y));
YSR(find(z<=20),:) = 1.5;
St = st.*ones(size(y));
% qc = qc.*reduction_f;
% sands
sgrc(S_ind,:) = ((0.029*qc(S_ind).*(sgp(S_ind)/101.3).^0.13)*ones(size(L))).*y(S_ind,:).^-0.38;
eta = qc(S_ind).*(101.3.*sgp(S_ind)).^-0.5;eta(isnan(eta)) = 0;
G = qc(S_ind)./(0.0203+0.00125*eta-1.213*10^-6.*eta.^2);
sgrd(S_ind,:) = (2.*2*10^-5/R.*G)*ones(size(L));
sgrf = sgrc(S_ind,:)+sgrd(S_ind,:);
tau(S_ind,:) = sgrf.*(tand(delta_cv(S_ind))*ones(size(L)));
ten(S_ind,:) = 0.9*(0.8*sgrc(S_ind,:)+sgrd(S_ind,:)).*(tand(delta_cv(S_ind))*ones(size(L)));
% silts or clay
Ivy = log10(St);
Kc = y(C_ind,:).^-0.2.*YSR(C_ind,:).^0.42.*(2.2+0.016.*YSR(C_ind,:)-0.87.*Ivy(C_ind,:));
sgrc(C_ind,:) = Kc.*(sgp(C_ind)*ones(size(L)));
tau(C_ind,:) = 0.8*sgrc(C_ind,:).*(tand(delta_cv(C_ind))*ones(size(L)));
ten(C_ind,:) = 0.8*sgrc(C_ind,:).*(tand(delta_cv(C_ind))*ones(size(L)));
% bearing capacity
[~,ind] = min(abs(z*ones(size(L))-ones(size(z))*L));
qc_D = movmean(qc',K)';
Qp = (Ic(ind)<2.6).*(qc_D(ind)*Ar*pi*(D.^2)/4) + (Ic(ind)>=2.6).*(qc(ind)*Ar*pi*(D.^2)/4);
temp_qc = qc*Ar;temp_qc(S_ind) = qc_D(S_ind)*Ar;
tau = tau.*cdg;
ten = ten.*cdg;
for i = 1:length(ind)
    Qf(i,:) = sum(pi*D*dz*tau(1:ind(i),i));
    Qt(i,:) = sum(pi*D*dz*ten(1:ind(i),i));
    Qz(1:ind(i),i) = cumsum(pi*D*dz*tau(1:ind(i),i)) + temp_qc(1:ind(i))*pi*(D.^2)/4;
    Qzt(1:ind(i),i) = cumsum(pi*D*dz*ten(1:ind(i),i));
    qbz(1:ind(i),i) =  temp_qc(1:ind(i));
    S_prc(i,:) = sum(S_ind<=ind(i))./ind(i);
    S_fc_prc(i,:) = sum(pi*D*dz*(S_ind<=ind(i)).*tau(S_ind,i))/Qf(i,:);
    S_ft_prc(i,:) = sum(pi*D*dz*(S_ind<=ind(i)).*ten(S_ind,i))/Qt(i,:);
end
Q = Qf+Qp;
end

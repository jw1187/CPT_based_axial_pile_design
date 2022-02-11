function [Q,Qf,Qp,Qt,Qz,Qzt,tau,qbz,S_prc,S_fc_prc,S_ft_prc]=UWA_full_liq(z,qc,Ic,sgp,D,WT,L,De)

dz = z(2)-z(1);
K = round(3*D/dz);
R = sqrt(D.^2-(D-2*WT)^2)/2;
cdg = 0.96;     % cyclic degradation factor
Ar = 1-((D-2*WT)^2./D^2);
S_ind = find(Ic<2.6);
C_ind = find(Ic>=2.6);
y = ((ones(size(z))*L-z*ones(size(L)))./D);
y(C_ind,:) = y(C_ind,:)*D/R;
y(find(y<=1)) = 1;
y(find(y<=2 & y>1)) = 2; 
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
% reduction_f = zeros(size(Ic));
delta_cv = (Ic<2.05).*30 + (Ic>=2.05 & Ic<2.6).*27 + (Ic>=2.6 & Ic<2.95).*25 + (Ic>=2.95).*20;
% sands
sgrc(S_ind,:) = (0.03*qc(S_ind).*Ar^0.3)*ones(size(L)).*y(S_ind,:).^-0.5;
tau(S_ind,:) = sgrc(S_ind,:).*(tand(delta_cv(S_ind)*ones(size(L))));
ten(S_ind,:) = 0.75*tau(S_ind,:);
% silts or clay
tau(C_ind,:) = (0.23*qc(C_ind).*(qc(C_ind)./sgp(C_ind)).^-0.15)*ones(size(L)).*y(C_ind,:).^-0.2;
tau(isnan(tau)) = 0;
tau(C_ind,:) = tau(C_ind,:).*(tand(delta_cv(C_ind)*ones(size(L))));
ten(C_ind,:) = tau(C_ind,:);
% bearing capacity
[~,ind] = min(abs(z*ones(size(L))-ones(size(z))*L));
qc_D = movmean(qc',K)';
Qp = (0.15+0.45*Ar).*qc_D(ind)*pi*(D.^2)/4;
tau = tau.*cdg.*(De*ones(size(L)));
ten = ten.*cdg.*(De*ones(size(L)));
for i = 1:length(ind)
    Qf(i,:) = sum(pi*D*dz*tau(1:ind(i),i));
    Qt(i,:) = sum(pi*D*dz*ten(1:ind(i),i));
    Qz(1:ind(i),i) = cumsum(pi*D*dz*tau(1:ind(i),i)) + qc_D(1:ind(i)).*(0.15+0.45*Ar).*pi*(D.^2)/4;
    Qzt(1:ind(i),i) = cumsum(pi*D*dz*ten(1:ind(i),i));
    qbz(1:ind(i),i) =  qc_D(1:ind(i)).*(0.15+0.45*Ar);
    S_prc(i,:) = sum(S_ind<=ind(i))./ind(i);
    S_fc_prc(i,:) = sum(pi*D*dz*(S_ind<=ind(i)).*tau(S_ind,i))/Qf(i,:);
    S_ft_prc(i,:) = sum(pi*D*dz*(S_ind<=ind(i)).*ten(S_ind,i))/Qt(i,:);
end
Q = Qf+Qp;
end
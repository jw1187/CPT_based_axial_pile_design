function [Q,Qf,Qp,Qt,Qz,Qzt,tau,qbz,S_prc,S_fc_prc,S_ft_prc]=Fugro_full(z,qc,Ic,sgp,D,WT,L,sg)

dz = z(2)-z(1);
K = round(3*D/dz);
R = sqrt(D.^2-(D-2*WT)^2)/2;
cdg = 0.96;     % cyclic degradation factor
y = ((ones(size(z))*L-z*ones(size(L)))./R);
y = (y.^-0.9).*(y>=4)+(y<4).*((4^-0.9)*y/4);y(find(isnan(y))) = 0;  %for compression
y1 = ((ones(size(z))*L-z*ones(size(L)))./R);    % for tension
y1(find(y1<4)) = 4;
y2 = ((ones(size(z))*L-z*ones(size(L)))); y2(find(y2<1e-10)) = 0;   % for clay

Ar = 1-((D-2*WT)^2./D^2);
S_ind = find(Ic<2.6);
C_ind = find(Ic>=2.6);
ks = zeros(size(y));
tau = zeros(size(y));
ten = zeros(size(y));
Qz = zeros(size(y));
Qzt = zeros(size(y));
qbz = zeros(size(y));
S_prc = zeros(size(L'));
S_fc_prc = zeros(size(L'));
S_ft_prc = zeros(size(L'));
% sands
tau(S_ind,:) = (0.08*qc(S_ind).*(sgp(S_ind)/101.3).^0.05)*ones(size(L)).*y(S_ind,:);tau(isnan(tau(S_ind,:))) = 0;
ten(S_ind,:) = (0.045*qc(S_ind).*(sgp(S_ind)/101.3).^0.15)*ones(size(L)).*y1(S_ind,:).^-0.85;ten(isnan(ten(S_ind,:))) = 0;
% silts or clay
temp_Qt = (qc(C_ind)-sg(C_ind))./sgp(C_ind); temp_Qt(isnan(temp_Qt) | isinf(temp_Qt)) = 0;temp_Qt(find(temp_Qt<0)) = 0;
ks(C_ind,:) = (0.16*(temp_Qt.^-0.4)*ones(size(L))).*y2(C_ind,:).^-0.3;ks(isinf(ks)) = 0;
ks(C_ind,:) = (ks(C_ind,:)>0.08).*0.08 + (ks(C_ind,:)<=0.08).*ks(C_ind,:);
tau(C_ind,:) = ks(C_ind,:).*(temp_Qt.*sgp(C_ind)*ones(size(L)));
ten(C_ind,:) = tau(C_ind,:);
tau = tau.*cdg;
ten = ten.*cdg;
% bearing capacity
[~,ind] = min(abs(z*ones(size(L))-ones(size(z))*L));
qc_D = movmean(qc',K)';
sg_D = movmean(sg',K)';
temp_qc = 8.5.*(qc_D/101.3).^0.5.*(Ar^0.25)*101.3;
temp_qc(C_ind) = 0.7*(qc_D(C_ind)-sg_D(C_ind));
Qp = temp_qc(ind)*pi*(D.^2)/4;
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


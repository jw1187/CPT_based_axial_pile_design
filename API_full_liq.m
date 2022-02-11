function [Q,Qf,Qp,Qt,Qz,Qzt,tau,qb]=API_full_liq(z,qc,sg,sgp,Ic,D,WT,L,De)
    
dz = z(2)-z(1);
% beta_b = [0.21,0.29,0.37,0.46,0.56];
R = sqrt(D.^2-(D-2*WT).^2)/2;
cdg = 0.96;     % cyclic degradation factor
Ar = 1-((D-2*WT)^2./D^2);
sup_f = [48,67,81,96,115];
Nq_b = [8,12,20,40,50];
sup_b = [2000,3000,5000,10000,12000];
Dr_bound  = [15,35,65,85];
S_ind = find(Ic<2.6);
C_ind = find(Ic>=2.6);
Qz = zeros(length(z),length(L));
Qzt = zeros(length(z),length(L));
beta = zeros(size(z));
Nq = zeros(size(z));
tau =  zeros(size(z));
 qb = zeros(size(z));
% beta_boundary
% sands
Dr = -98+66*log10(qc(S_ind)/9.8./(sgp(S_ind)/9.8).^0.5);Dr(isnan(Dr)) = 0;
ind1 = find(Dr<Dr_bound(1)); beta(S_ind(ind1)) = 0.8*tand(15);Nq(S_ind(ind1)) = Nq_b(1);
ind2 = find(Dr<Dr_bound(2) & Dr>=Dr_bound(1)); beta(S_ind(ind2)) = 0.8*tand(20);Nq(S_ind(ind2)) = Nq_b(2);
ind3 = find(Dr<Dr_bound(3) & Dr>=Dr_bound(2)); beta(S_ind(ind3)) = 0.8*tand(25);Nq(S_ind(ind3)) = Nq_b(3);
ind4 = find(Dr<Dr_bound(4) & Dr>=Dr_bound(3)); beta(S_ind(ind4)) = 0.8*tand(30);Nq(S_ind(ind4)) = Nq_b(4);
ind5 = find(Dr>=Dr_bound(4)); beta(S_ind(ind5)) = 0.8*tand(35);Nq(S_ind(ind5)) =Nq_b(5);
%clay
su = (qc(C_ind)-sg(C_ind))/15; su(find(su<0)) = 0;%su±j­¢=0
ind6 = find(su./sgp(C_ind)>1);beta(C_ind(ind6)) = 0.5*(su(ind6)./sgp(C_ind(ind6))).^-0.25;beta(C_ind(ind6)) = ~isinf(beta(C_ind(ind6))).*beta(C_ind(ind6));
ind7 = find(su./sgp(C_ind)<=1);beta(C_ind(ind7)) = 0.5*(su(ind7)./sgp(C_ind(ind7))).^-0.5;beta(C_ind(ind7)) = ~isinf(beta(C_ind(ind7))).*beta(C_ind(ind7));
ind8 = find(su<=0 | sgp(C_ind)<=0); beta(C_ind(ind8)) = 0;
qb(C_ind) =9*su;
%tau
tau(S_ind) = beta(S_ind).*sgp(S_ind);
tau(S_ind(ind1)) = (tau(S_ind(ind1))>sup_f(1)).*sup_f(1)+(tau(S_ind(ind1))<=sup_f(1)).*tau(S_ind(ind1));
tau(S_ind(ind2)) = (tau(S_ind(ind2))>sup_f(2)).*sup_f(2)+(tau(S_ind(ind2))<=sup_f(2)).*tau(S_ind(ind2));
tau(S_ind(ind3)) = (tau(S_ind(ind3))>sup_f(3)).*sup_f(3)+(tau(S_ind(ind3))<=sup_f(3)).*tau(S_ind(ind3));
tau(S_ind(ind4)) = (tau(S_ind(ind4))>sup_f(4)).*sup_f(4)+(tau(S_ind(ind4))<=sup_f(4)).*tau(S_ind(ind4));
tau(S_ind(ind5)) = (tau(S_ind(ind5))>sup_f(5)).*sup_f(5)+(tau(S_ind(ind5))<=sup_f(5)).*tau(S_ind(ind5));
tau(C_ind) = beta(C_ind).*su;
tau = tau.*cdg.*De;
%qb
qb(S_ind) = Nq(S_ind).*sgp(S_ind);
qb(S_ind(ind1)) = (qb(S_ind(ind1))>sup_b(1)).*sup_b(1)+(qb(S_ind(ind1))<=sup_b(1)).*qb(S_ind(ind1));
qb(S_ind(ind2)) = (qb(S_ind(ind2))>sup_b(2)).*sup_b(2)+(qb(S_ind(ind2))<=sup_b(2)).*qb(S_ind(ind2));
qb(S_ind(ind3)) = (qb(S_ind(ind3))>sup_b(3)).*sup_b(3)+(qb(S_ind(ind3))<=sup_b(3)).*qb(S_ind(ind3));
qb(S_ind(ind4)) = (qb(S_ind(ind4))>sup_b(4)).*sup_b(4)+(qb(S_ind(ind4))<=sup_b(4)).*qb(S_ind(ind4));
qb(S_ind(ind5)) = (qb(S_ind(ind5))>sup_b(5)).*sup_b(5)+(qb(S_ind(ind5))<=sup_b(5)).*qb(S_ind(ind5));
%
[~,ind] = min(abs(z*ones(size(L))-ones(size(z))*L));
for i = 1:length(ind)
    Qf(i,:) = sum(pi*D*dz*tau(1:ind(i)));
    Qt(i,:) = sum(pi*D*dz*tau(1:ind(i)));
    Qp(i,:) = min(qb(ind(i))*Ar*pi*(D.^2)/4 + sum(pi*(D-2*WT)*dz*tau(1:ind(i))), qb(ind(i))*pi*(D.^2)/4) ;
    Qz(1:ind(i),i) = cumsum(pi*D*dz*tau(1:ind(i))) + min([qb(1:ind(i)).*pi*(D.^2)/4,qb(1:ind(i))*Ar*pi*(D.^2)/4 + cumsum(pi*(D-2*WT)*dz*tau(1:ind(i)))]')';
    Qzt(1:ind(i),i) = cumsum(pi*D*dz*tau(1:ind(i)));
%     S_prc(i,:) = sum(S_ind<=ind(i))./ind(i);
%     S_fc_prc(i,:) = sum(pi*D*dz*(S_ind<=ind(i)).*tau(S_ind,i))/Qf(i,:);
%     S_ft_prc(i,:) = sum(pi*D*dz*(S_ind<=ind(i)).*ten(S_ind,i))/Qt(i,:);
end
Q = Qf+Qp;


end

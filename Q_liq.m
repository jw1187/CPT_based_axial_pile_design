function [Q,Qt] = Q_liq(d,qt,Ic,sg,sgp,D,WT,L,st,De)
for nnn = 1:3
    [Q(:,1,nnn),~,~,Qt(:,1,nnn)]=ICP05_full_liq(d,qt,Ic,sgp,D,WT,L,st,De(:,nnn));
    [Q(:,2,nnn),~,~,Qt(:,2,nnn)]=UWA_full_liq(d,qt,Ic,sgp,D,WT,L,De(:,nnn));
    [Q(:,3,nnn),~,~,Qt(:,3,nnn)]=NGI_full_liq(d,qt,Ic,sg,sgp,D,WT,L,De(:,nnn));
    [Q(:,4,nnn),~,~,Qt(:,4,nnn)]=Fugro_full_liq(d,qt,Ic,sgp,D,WT,L,sg,De(:,nnn));
    [Q(:,5,nnn),~,~,Qt(:,5,nnn)]=API_full_liq(d,qt,sg,sgp,Ic,D,WT,L,De(:,nnn));
end
end
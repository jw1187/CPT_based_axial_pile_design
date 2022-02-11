function beta = BetaC(sai)
    beta = zeros(size(sai));
    ind6 = find(sai>1);beta((ind6)) = 0.5*(sai(ind6)).^-0.25;beta((ind6)) = ~isinf(beta((ind6))).*beta((ind6));
    ind7 = find(sai<=1);beta((ind7)) = 0.5*(sai(ind7)).^-0.5;beta((ind7)) = ~isinf(beta((ind7))).*beta((ind7));
    beta(find(sai==0)) = 0;
end
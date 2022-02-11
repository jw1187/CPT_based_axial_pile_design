function beta = BetaS(Dr)
    Dr_bound  = [15,35,65,85];
    beta = zeros(size(Dr));
    ind1 = find(Dr<Dr_bound(1)); beta(ind1) = 0.8*tand(15);
    ind2 = find(Dr<Dr_bound(2) & Dr>=Dr_bound(1)); beta((ind2)) = 0.8*tand(20);
    ind3 = find(Dr<Dr_bound(3) & Dr>=Dr_bound(2)); beta((ind3)) = 0.8*tand(25);
    ind4 = find(Dr<Dr_bound(4) & Dr>=Dr_bound(3)); beta((ind4)) = 0.8*tand(30);
    ind5 = find(Dr>=Dr_bound(4)); beta((ind5)) = 0.8*tand(35);
end
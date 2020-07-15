function res = BCs_eta_left_gamma_right_incom(yl,yr)
global deltatilde
res = zeros(14,1);
res(1) = yl(2);
res(2) = yl(4);
res(3) = yr(11);
res(4) = yr(13);
res(5) = yr(1)-yl(6);
res(6) = yr(2)-yl(7);
res(7) = yr(3)-yl(8);
res(8) = yr(4)-yl(9);
res(9) = yr(6) - yl(10);
res(10) = yr(7) - yl(11);
res(11) = yr(8) - yl(12);
res(12) = yr(9) - yl(13);
res(13) = yr(5) + deltatilde/yl(6)/yl(8);
res(14) = yl(14) + deltatilde/yr(6)/yr(8);
end
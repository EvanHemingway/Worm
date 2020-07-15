function res = BCs_incom_gam(yl,yr)
global deltatilde
res = zeros(9,1);
res(1) = yl(2);
res(2) = yl(4);
res(3) = yr(6);
res(4) = yr(8);
res(5) = yr(1)-yl(5);
res(6) = yr(2)-yl(6);
res(7) = yr(3)-yl(7);
res(8) = yr(4)-yl(8);
res(9) = yl(9) + deltatilde/yr(1)/yr(3);
end
function res = BCs_incom_eta(yl,yr)
global deltatilde
res = zeros(9,1);
res(1) = yl(2);
res(2) = yl(4);
res(3) = yr(7);
res(4) = yr(9);
res(5) = yr(1)-yl(6);
res(6) = yr(2)-yl(7);
res(7) = yr(3)-yl(8);
res(8) = yr(4)-yl(9);
res(9) = yr(5) + deltatilde/yl(6)/yl(8);
end
function res = BCs_long_con_incom(yl,yr)
res = zeros(4,1);
res(1) = yl(2);
res(2) = yl(4);
res(3) = yr(2);
res(4) = yr(4);
end
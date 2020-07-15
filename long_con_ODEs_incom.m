function yp = resting_ODEs_incom(zeta,y)
global itilde mtilde
yp = zeros(4,1);
yp(1) = y(2);
yp(2) = mtilde*y(1)*y(3)-deltatilde/y(1)^2/y(3)+itilde*(2*y(1)+y(3)-3);
yp(3) = y(4);
yp(4) = -deltatilde/y(1)/y(3)^2+itilde*(2*y(3)+y(1)-3);
end
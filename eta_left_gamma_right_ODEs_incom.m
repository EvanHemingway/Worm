function yp = eta_left_gamma_right_ODEs_incom(zeta,y)
global itilde deltatilde ctilde mtilde etatilde gamtilde
yp = zeros(14,1);
yp(1) = etatilde*y(2);
yp(2) = etatilde*(ctilde + mtilde*y(1)*y(3)+y(5)/y(1)+itilde*(2*(y(1)-1)+(y(3)-1)));
yp(3) = etatilde*y(4);
yp(4) = etatilde*(ctilde+ y(5)/y(3)+itilde*(2*(y(3)-1)+(y(1)-1)));
yp(5) = - etatilde^2*y(5)*(y(2)/y(1)+y(4)/y(3));
yp(6) = (gamtilde-etatilde)*y(7);
yp(7) = (gamtilde-etatilde)*(mtilde*y(6)*y(8)-deltatilde/y(6)^2/y(8)+itilde*(2*(y(6)-1)+(y(8)-1)));
yp(8) = (gamtilde-etatilde)*y(9);
yp(9) = (gamtilde-etatilde)*(-deltatilde/y(6)/y(8)^2+itilde*(2*(y(8)-1)+(y(6)-1)));
yp(10) = (1-gamtilde)*y(11);
yp(11) = (1-gamtilde)*(ctilde + mtilde*y(10)*y(12)+y(14)/y(10)+itilde*(2*(y(10)-1)+(y(12)-1)));
yp(12) = (1-gamtilde)*y(13);
yp(13) = (1-gamtilde)*(ctilde+ y(14)/y(12)+itilde*(2*(y(12)-1)+(y(10)-1)));
yp(14) = - (1-gamtilde)^2*y(14)*(y(11)/y(10)+y(13)/y(12));
end
function yp = just_eta_ODEs_incom(zeta,y)
global itilde deltatilde ctilde mtilde etatilde gamtilde
yp = zeros(9,1);
yp(1) = etatilde*y(2);
yp(2) = etatilde*(ctilde + mtilde*y(1)*y(3)+y(5)/y(1)+itilde*(2*(y(1)-1)+(y(3)-1)));
yp(3) = etatilde*y(4);
yp(4) = etatilde*(ctilde+ y(5)/y(3)+itilde*(2*(y(3)-1)+(y(1)-1)));
yp(5) = - etatilde^2*y(5)*(y(2)/y(1)+y(4)/y(3));
yp(6) = (1-etatilde)*y(7);
yp(7) = (1-etatilde)*(mtilde*y(6)*y(8)-deltatilde/y(6)^2/y(8)+itilde*(2*(y(6)-1)+(y(8)-1)));
yp(8) = (1-etatilde)*y(9);
yp(9) = (1-etatilde)*(-deltatilde/y(6)/y(8)^2+itilde*(2*(y(8)-1)+(y(6)-1)));
end
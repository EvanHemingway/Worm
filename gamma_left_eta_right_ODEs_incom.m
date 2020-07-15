function yp = gamma_left_eta_right_ODEs_incom(zeta,y)
global itilde deltatilde ctilde mtilde etatilde gamtilde
yp = zeros(13,1);
yp(1) = gamtilde*y(2);
yp(2) = gamtilde*(mtilde*y(1)*y(3)-deltatilde/y(1)^2/y(3)+itilde*(2*(y(1)-1)+(y(3)-1)));
yp(3) = gamtilde*y(4);
yp(4) = gamtilde*(-deltatilde/y(1)/y(3)^2+itilde*(2*(y(3)-1)+(y(1)-1)));
yp(5) = (etatilde-gamtilde)*y(6);
yp(6) = (etatilde-gamtilde)*(ctilde + mtilde*y(5)*y(7)+y(9)/y(5)+itilde*(2*(y(5)-1)+(y(7)-1)));
yp(7) = (etatilde-gamtilde)*y(8);
yp(8) = (etatilde-gamtilde)*(ctilde+ y(9)/y(7)+itilde*(2*(y(7)-1)+(y(5)-1)));
yp(9) = - (etatilde-gamtilde)^2*y(9)*(y(6)/y(5)+y(8)/y(7));
yp(10) = (1-etatilde)*y(11);
yp(11) = (1-etatilde)*(mtilde*y(10)*y(12)-deltatilde/y(10)^2/y(12)+itilde*(2*(y(10)-1)+(y(12)-1)));
yp(12) = (1-etatilde)*y(13);
yp(13) = (1-etatilde)*(-deltatilde/y(10)/y(12)^2+itilde*(2*(y(12)-1)+(y(10)-1)));
end
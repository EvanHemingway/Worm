function yp = just_gamma_ODEs_incom(zeta,y)
global itilde deltatilde ctilde mtilde etatilde gamtilde
yp = zeros(9,1);
yp(1) = gamtilde*y(2);
yp(2) = gamtilde*(mtilde*y(1)*y(3)-deltatilde/y(1)^2/y(3)+itilde*(2*(y(1)-1)+(y(3)-1)));
yp(3) = gamtilde*y(4);
yp(4) = gamtilde*(-deltatilde/y(1)/y(3)^2+itilde*(2*(y(3)-1)+(y(1)-1)));
yp(5) = (1-gamtilde)*y(6);
yp(6) = (1-gamtilde)*(ctilde + mtilde*y(5)*y(7)+y(9)/y(5)+itilde*(2*(y(5)-1)+(y(7)-1)));
yp(7) = (1-gamtilde)*y(8);
yp(8) = (1-gamtilde)*(ctilde+y(9)/y(7)+itilde*(2*(y(7)-1)+(y(5)-1)));
yp(9) = - (1-gamtilde)^2*y(9)*(y(6)/y(5)+y(8)/y(7));
end
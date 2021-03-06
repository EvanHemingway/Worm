function [mu1,mu1p,mu2,mu2p,xi]=long_con_incom()
global ell N

%Allocate a uniform distribution of zeta points

zetamesh = linspace(0,1,N);

solinit = bvpinit(zetamesh,@IG_long_con_incom);

%options = bvpset('Stats','on');
sol = bvp4c(@long_con_ODEs_incom,@BCs_long_con_incom,solinit);
%sol.stats

z = sol.x;
mu1 = sol.y(1,:);
mu1p = sol.y(2,:)/ell;
mu2 = sol.y(3,:);
mu2p = sol.y(4,:)/ell;
xi = z*ell;
end

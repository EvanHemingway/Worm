function [lam,mu1,mu1p,mu2,mu2p,xi]=just_eta_incom()
global ell N
global etatilde Dbar kG

%Allocate a uniform distribution of zeta points for left and right domains

zetamesh = linspace(0,1,N);

solinit = bvpinit(zetamesh,@IG_incom_eta);

%options = bvpset('Stats','on');
sol = bvp4c(@just_eta_ODEs_incom,@BCs_incom_eta,solinit);
%sol.stats

z = sol.x;
zl = z*etatilde;
zr = z*(1-etatilde) + etatilde;
mu1l = sol.y(1,:);
mu1pl = sol.y(2,:)/ell;
mu2l = sol.y(3,:);
mu2pl = sol.y(4,:)/ell;
laml = kG*sol.y(5,:)/ell^2;
mu1r = sol.y(6,:);
mu1pr = sol.y(7,:)/ell;
mu2r = sol.y(8,:);
mu2pr = sol.y(9,:)/ell;
lamr = -Dbar./(mu1r.*mu2r);

z = [zl,zr(2:end)];
xi = z*ell;
mu1 = [mu1l, mu1r(2:end)];
mu2 = [mu2l, mu2r(2:end)];
mu1p = [mu1pl, mu1pr(2:end)];
mu2p = [mu2pl, mu2pr(2:end)];
lam = [laml, lamr(2:end)];
end
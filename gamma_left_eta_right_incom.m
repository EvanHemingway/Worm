function [lam,mu1,mu1p,mu2,mu2p,xi]=gamma_left_eta_right_incom()
global ell N
global etatilde gamtilde kG Dbar

zetamesh = linspace(0,1,N);

solinit = bvpinit(zetamesh,@IG_nonempty_incom);

%options = bvpset('Stats','on');
sol = bvp4c(@gamma_left_eta_right_ODEs_incom,@BCs_nonempty_incom,solinit);
%sol.stats

z = sol.x;
zl = z*gamtilde;
zm = z*(etatilde-gamtilde)+gamtilde;
zr = z*(1-etatilde) + etatilde;
mu1l = sol.y(1,:);
mu1pl = sol.y(2,:)/ell;
mu2l = sol.y(3,:);
mu2pl = sol.y(4,:)/ell;
laml = -Dbar./(mu1l.*mu2l);
mu1m = sol.y(5,:);
mu1pm = sol.y(6,:)/ell;
mu2m = sol.y(7,:);
mu2pm = sol.y(8,:)/ell;
lamm = kG*sol.y(9,:)/ell^2;
mu1r = sol.y(10,:);
mu1pr = sol.y(11,:)/ell;
mu2r = sol.y(12,:);
mu2pr = sol.y(13,:)/ell;
lamr = -Dbar./mu1r./mu2r;

z = [zl,zm(2:end),zr(2:end)];
xi = z*ell;
mu1 = [mu1l, mu1m(2:end), mu1r(2:end)];
mu2 = [mu2l, mu2m(2:end), mu2r(2:end)];
mu1p = [mu1pl, mu1pm(2:end), mu1pr(2:end)];
mu2p = [mu2pl, mu2pm(2:end), mu2pr(2:end)];
lam = [laml, lamm(2:end), lamr(2:end)];
end
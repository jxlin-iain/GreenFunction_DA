function G = tensorG(k0,r1x,r1y,r1z,r2x,r2y,r2z)
% Green function for two-dimensional periodic array
R12 = sqrt((r1x-r2x).^2+(r1y-r2y).^2+(r1z-r2z).^2);
er = [r2x-r1x,r2y-r1y,r2z-r1z]./R12;
G = ((1./R12 + 1i./(k0*R12.^2) - 1./(k0^2*R12.^3)) * [1,0,0;0,1,0;0,0,1] ...
    + (-1./R12 - 1i*3./(k0*R12.^2) + 3./(k0^2*R12.^3)) ...
* [er(1)*er(1),er(1)*er(2),er(1)*er(3); ...
er(2)*er(1),er(2)*er(2),er(2)*er(3);...
er(3)*er(1),er(3)*er(2),er(3)*er(3)]) .* exp(1i*k0*R12) / 4/ pi;

% exp(1i*k0*R12)




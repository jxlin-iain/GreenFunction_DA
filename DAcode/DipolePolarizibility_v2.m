% This code is to calculate the effective polarizability of single dipole
% and try to reproduce the macrostropic absorption. Make sure the result is
% consistent with CDM.
% version 2.0: eps = [1-e,1-e,1-2.9^2]
% extinction / absorption

clear all
clc;
lda = linspace(600,900,500)*1e-9;
omega = 299792458 *2 *pi ./ lda;
k0 = 2*pi./lda;
E0 = [sqrt(1.3707e17);0;0];
n = 5000;
spacing = 5*1e-9;
a = spacing;
L = 1 * 1E-9;
coe = 2/sqrt(3);
Gtot = zeros(3,3,length(lda));
Gcol = zeros(3,3,length(lda));
omega0 = 299792458*k0;
eV = 299792458 * 6.58211899e-16 * (2 * pi) ./ lda;
eps0 = 8.8541878188*1e-12;%8.8541878188*1e-12 F⋅m−1 vacuum dielectric constant
fA1 = 0.8; %1.9[eV]
omegaA1 = 1.6674; %2.014[eV] 
gammaA1 = 0.035;
epsinf = 1;
eps1 = epsinf + fA1./(omegaA1^2-eV.^2-1i*gammaA1.*eV);
% alphaE = eps0 * 4*pi*R^3*(epsd-1)./(epsd+2);
n1 = 1.15; % refractive index of the surrounding
kd = n1.*k0;
%% full
% sum
a = spacing;
c = 50;
G0l = zeros(3,3,length(lda));
G2inf = zeros(3,3,length(lda));
Gtot = zeros(3,3,length(lda));
NarrayCut = c;
element = generate_hexagonal_structure(floor(NarrayCut*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);
mask = R < NarrayCut;
 % mask = R < Narray;
x = x(mask);
y = y(mask);
element = a*[x(:),y(:)];
origin_idx = find(element(:,1) == 0 & element(:,2) == 0, 1);

if ~isempty(origin_idx)
    element = circshift(element, -(origin_idx - 1));
end
    total_elements = size(element, 1);  % Total iterations for the inner loop
for n = 1:length(lda)
    for ii = 2:total_elements
        Gij = tensorG(kd(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0l(:,:,n) = G0l(:,:,n) + Gij;
    end
    G2inf(1,1,n) = ((1/4).*a.^(-3).*c.^(-1).*kd(n).^(-2).*(exp(1) ...
  .^(sqrt(-1).*a.*c.*kd(n)).*(1+sqrt(-1).*a.*c.*kd(n))))*coe;
    G2inf(2,2,n) = G2inf(1,1,n);
    G2inf(3,3,n) = (1/2).*a.^(-3).*c.^(-1).*kd(n).^(-2).*(sqrt(-1).*exp(1).^( ...
  sqrt(-1).*a.*c.*kd(n)).*(sqrt(-1)+a.*c.*kd(n)))*coe;
    Gtot(:,:,n) = G0l(:,:,n) + G2inf(:,:,n);
end
%% coulomb
Gcol = zeros(3,3,length(lda));

for n = 1:length(lda)
    Gcolxx(n) = 1/4*a.^(-3).*k0(n).^(-2).*(exp(1i.*a.*k0(n))-1i.*a.*k0(n).*(ExpintEi(1i.*a.*k0(n))-1i*pi));
end
for n = 1:length(lda)
Gcol(:,:,n) = [Gcolxx(n),0,0;0,Gcolxx(n),0;0,0,-2*Gcolxx(n)]*coe;
end

%% single-dipole polarizability

V = spacing^2 *sqrt(3)/2 * L;
alphaeff = zeros(3,3,length(lda));
alphaE = zeros(3,3,length(lda));
sigma = zeros(length(lda),1);
for kk = 1:length(k0)
% alphaeff(:,:,kk) = vpa([eps1(kk) - 1,0,0;0,eps1(kk)-1,0;0,0,1-1/2.9^2]);
alphaeff(:,:,kk) = vpa([eps1(kk) - 1,0,0;0,eps1(kk)-1,0;0,0,2.9^2-1]);
alphaE(:,:,kk) = vpa(eye(3)/ (1 / V * eye(3)/alphaeff(:,:,kk) + k0(kk)^2 * Gcol(:,:,kk)));
% Pinf(:,kk) = eps0 * eye(3)/(eye(3)/alphaE(:,:,kk) - k0(kk)^2 * Gtot(:,:,kk)) * (transfermatrix(E0,lda(kk),1.15^2,eps1(kk),2.25,1e-9) * [1;0;0]);
Pinf(:,kk) = eps0 * eye(3)/(eye(3)/alphaE(:,:,kk) - k0(kk)^2 * Gtot(:,:,kk)) * E0;
sigma(kk) = - omega0(kk) * sum(imag(conj(Pinf(:,kk)).*eye(3)/(alphaE(:,:,kk))/eps0 * Pinf(:,kk)));
end

%% transfer matrix for Reflectivity & Fersnel reflection
theta0 = 0; % incident angle, unit: degree
ntotal = 1000;
 %%%%%% only for the normal incident direction
z0 = linspace(0,1,ntotal)*1e-9;
E1 = zeros(3,length(z0),length(lda));
R = zeros(1,length(lda));
for jj = 1:length(lda)

k0 = 2*pi/lda(jj);% wavevector in vacuum

% permittivity of layer1  ; e.g., MoS2

eps2 = 1.5;% permittivity of substrate  ; e.g., Si
epsd = 1.15; %1.15
n0 = sqrt(epsd); % refractive index of environment ; e.g., air
n1 = sqrt(eps1(jj)); 
n2 = sqrt(eps2);

theta1 = asind(n0./n1 .* sind(theta0));
theta2 = asind(n1./n2 .* sind(theta1));
k0z = n0.*k0.*cosd(theta0);
k1z = n1.*k0.*cosd(theta1);
k2z = n2.*k0.*cosd(theta2);
kappa1 = k0z./k1z;
kappa2 = k1z./k2z;
eta1 = eps1(jj)/epsd;
eta2 = eps2./eps1(jj);

Phi1 = [exp(-1i*k1z*L),0;0,exp(1i*k1z*L)]; % propagation matrix
T01 = 0.5*[1+kappa1.*eta1,1-kappa1.*eta1;1-kappa1.*eta1,1+kappa1.*eta1]; % transmission matrix
T02 = 0.5*[1+kappa2.*eta2,1-kappa2.*eta2;1-kappa2.*eta2,1+kappa2.*eta2];

E2i = 0.06594/0.14273;% amplitude of the incident field in last medium
Eind = T01 * Phi1 * T02 * [E2i;0]; 

R(jj) = abs(Eind(2))^2./abs(Eind(1))^2;

kx = n0 * k0.*sind(theta0);
t01(jj) = 2 * k0z * sqrt(epsd*eps1(jj)) / (eps1(jj)*k0z+epsd*k1z);
t10(jj) = 2 * k1z * sqrt(epsd*eps1(jj)) / (eps1(jj)*k0z+epsd*k1z);

r10(jj) = (epsd*k1z - eps1(jj) * k0z)/(epsd*k1z + eps1(jj) * k0z);
r12(jj) = (eps2*k1z - eps1(jj) * k2z)/(eps2*k1z + eps1(jj) * k2z);
r01(jj) = -r10(jj);
r21(jj) = -r12(jj);
uin0 = [k0z/sqrt(k0z^2+kx^2);kx/sqrt(k0z^2+kx^2)];
uin1 = [k1z/sqrt(k1z^2+kx^2);kx/sqrt(k1z^2+kx^2)];
ur = [-k0z/sqrt(k0z^2+kx^2);kx/sqrt(k0z^2+kx^2)];
z = 0.5 * 1e-9;
x = 0;
% E1 = norm(E0)*t01*exp(1i*kx*x-1i*k1z*(z-L)).*uin0/(1-r12(jj)*r10(jj)*exp(1i*2*k1z*L)) ...
% +norm(E0)*t01*exp(1i*kx*x+1i*k1z*z)* r12(jj)*exp(1i*k1z*L)*uin1/(1-r12(jj)*r10(jj)*exp(2*1i*k1z*L));

 E0r(1:2,jj) = norm(E0)*exp(1i*kx*x+1i*k0z*(z-L))*...
    (r01(jj)+t01(jj)*t10(jj)*r12(jj)*exp(1i*2*k1z*L)/(1-r12(jj)*r10(jj)*exp(1i*2*k1z*L)))*ur;

 E0i(1:2,jj) =  norm(E0)*exp(1i*kx*x-1i*k0z*(z-L))*uin0;

  E1(1,:,jj) = norm(E0)*t01(jj).*exp(1i*kx*x-1i*k1z*(z0-L))./(1-r12(jj)*r10(jj).*exp(1i*2*k1z*L)) ...
+norm(E0)*t01(jj).*exp(1i*kx*x+1i*k1z*z0).* r12(jj).*exp(1i*k1z*L)./(1-r12(jj)*r10(jj)*exp(2*1i*k1z*L));

S = (800 * 1e-9)^2;
Vint = S*1e-9;

Pabs(jj) = real(-1i*omega(jj)*sum((eps0*[eps1(jj),0,0;0,eps1(jj),0;0,0,2.9^2]*(conj(E1(:,jj)).*(E1(:,jj))))))/ntotal * Vint; % W
%ewfd.Jx = ewfd.iomega*(epsilon0_const*ewfd.Ex+epsilon0_const*(ewfd.epsilonrxx*ewfd.Ex-ewfd.Ex))
end
Rfn = (abs(E0r(1,:)).^2)./abs(E0i(1,:)).^2;
% figure(1)
% yyaxis left
% plot(eV,R)
% yyaxis right
% plot(eV,Rfn,'--')

%%

figure(1)
yyaxis left
plot(eV,squeeze(imag(alphaE(1,1,:))))
ylabel('Im(\alpha_{dip})')
yyaxis right
plot(eV,squeeze(imag(alphaeff(1,1,:))))
ylabel('Im(\alpha_{mat}')

figure(2)
yyaxis left
plot(eV,sigma)
ylabel('Infinite absorption')
yyaxis right
plot(eV,Pabs)
ylabel('layer absorption')

figure(3)
yyaxis left
plot(eV,sigma)
ylabel('infinite absorpion')
yyaxis right
plot(eV,squeeze(imag(alphaeff(1,1,:))))
ylabel('Im(\alpha_{mat}')

figure(4)
yyaxis left
plot(eV,Pabs/max(Pabs))
ylabel('normalized layer absorption')
yyaxis right
plot(eV,R/max(R))
ylabel('normalzied layer reflectivity')

Imalphaeffxx = double(squeeze(imag(alphaeff(1,1,:))));
ImalphaDipxx = double(squeeze(imag(alphaE(1,1,:))));
ext = double(sigma);
% result = real([eV',ext,Imalphaeffxx,ImalphaDipxx,Pabs']);
result = real([eV',ext,ImalphaDipxx]);
out = Pabs/max(Pabs);
% clc
% clear alphaeff;
% clear alphaE;
% clear alphaback;
% for kk = 1:length(k0)
% alphaeff(:,:,kk) = vpa(eps0 * [epsd(kk) - 1,0,0;0,epsd(kk)-1,0;0,0,1-1/2.9^2]);
% alphaE(:,:,kk) = vpa(eps0 * eye(3)  / (eps0 * eye(3) / alphaeff(:,:,kk) + k0(kk)^2 * G0tot(:,:,kk)));
% alphaback(:,:,kk) = vpa(eps0 * eye(3) / (eps0 *eye(3) / alphaE(:,:,kk) - k0(kk)^2 * G0tot(:,:,kk)));
% end
% figure(1)
% yyaxis left
% plot(eV,squeeze(imag(alphaback(1,1,:))))
% yyaxis right
% plot(eV,squeeze(imag(alphaeff(1,1,:))))
% title('ima(a)')
% figure(2)
% plot(eV,squeeze(imag(alphaE(1,1,:))))
% title('ima(aE)')

out = [eV;Pabs]';
out1 = real(squeeze(alphaE(1,1,:)));
out2 = imag(squeeze(alphaE(1,1,:)));
out = [sigma,out1,out2];

Ntot = Pabs*V./sigma';
figure(5)
plot(eV,Ntot)
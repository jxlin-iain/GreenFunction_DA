clear all
clc;
addpath('E:\ljx\Spain\rod_monolayer\ratio3');
import com.comsol.model.*
import com.comsol.model.util.*
model = mphload('WSe2_mesh_715nm.mph');
% Choose ur parameter

%% constant
eps0 = 8.8541878188*1e-12; %8.8541878188*1e-12 F⋅m−1 vacuum dielectric constant
econst = 1.602176634 * 1e-19; %C
hbar = 1.054571817*1e-34 ; %J`s
cconst = 299792458;

%% parameters
% dipoles properties
L = 1 * 1e-9; % height of 2D mat.
Narray = 20; % element number generator
spacing = 5*1e-9; %distance between dipoles
fA1 = 0.9; %1.9[eV]
omegaA1 = 1.6674; %2.014[eV]  1.6674 1.7072
gammaA1 = 0.035; % 0.0281 0.018
epsinf = 1; 
n1 = 1; % refractive index of air
Amplifier = 1; % amplitude for polarizability
shift =  0; %postion of the origin point of array
A = 1.79e6; % coefficient for Gmnp(R,R0) 1.79e6
flag = true;
configs = struct(...
    'use_R', {1, 2, 3, 4}, ...
    'param_name', {'spacing', 'Narray', 'Amplifier','Repeating'}, ...
    'paramArray', {[1,2,3,4,5]*1e-9, [20,30,40], linspace(10,100,10),[1:1:200]});
use_R = 4; 
lda = linspace(650,800,100) * 1e-9; % incident wavelength
% lda = linspace(400,700,100) * 1e-9;
E0free = [sqrt(1.3707e17);0;0];%V/m amplitude of incident electric field in free space

coe = 2/sqrt(3); % 1 for square, 2/sqrt(3) for Hexago
quarter = 1;
%%
k0 = 2*pi./lda;
kd = n1*k0;
omega0 = 299792458 * k0;
eV = 299792458 * 6.58211899e-16 * (2 * pi) ./ lda;
epsd = epsinf + fA1./(omegaA1^2-eV.^2-1i*gammaA1.*eV); % dielectric function for dipole

%% situation selection

current_config = configs([configs.use_R] == use_R);
paramArray = current_config.paramArray;
totalIterations = length(lda)*length(paramArray); % Total iterations for overall progress tracking
iterationCount = 0; % Counter for tracking overall progress
Qext_tot = zeros(length(paramArray),length(lda));
Qabs_tot = zeros(length(paramArray),length(lda));
Qext_dip = zeros(length(paramArray),length(lda));
Qext_comp = zeros(length(paramArray),length(lda));
                    ImGmnpxx = zeros(length(lda),1);
                    ImGijxx = zeros(length(lda),1);
                    ReGmnpxx = zeros(length(lda),1);
                    ReGijxx = zeros(length(lda),1);

for xx = 1:length(paramArray)
% Qext_dip(xx,1)
eval([current_config.param_name '=' num2str(paramArray(xx)) ';']);
flag = true;

if flag
%% Infinite summation of full dipole-dipole green function
a = spacing;
c = 50; % start integraling radius default 50
G0l = zeros(3,3,length(lda));
G2inf = zeros(3,3,length(lda));
Gtot = zeros(3,3,length(lda));
NarrayCut = c;
if coe==2/sqrt(3)
element = generate_hexagonal_structure(floor(NarrayCut*2/sqrt(3)));
elseif coe==1
    [x, y] = meshgrid(-(NarrayCut):1:(NarrayCut), ...
                  -(NarrayCut):1:(NarrayCut));
element = [x(:),y(:)];
end
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);
mask = R < NarrayCut;
 % mask = R < Narray;
x = x(mask);
y = y(mask);
element = a*[x(:),y(:)];
origin_idx = find(element(:,1) == 0 & element(:,2) == 0, 1);
% scatter(element(:,1),element(:,2))
if ~isempty(origin_idx)
    element = circshift(element, -(origin_idx - 1));
end
    total_elements = size(element, 1);  % Total iterations for the inner loop
for n = 1:length(lda)
G0l_ini = tensorG_vector(kd(n),zeros(total_elements,1),zeros(total_elements,1),zeros(total_elements,1),element(:,1),element(:,2),zeros(total_elements,1));
G0l(:,:,n) = sum(G0l_ini,3);
    G2inf(1,1,n) = ((1/4).*a.^(-3).*c.^(-1).*kd(n).^(-2).*(exp(1) ...
  .^(sqrt(-1).*a.*c.*kd(n)).*(1+sqrt(-1).*a.*c.*kd(n))))*coe;
    G2inf(2,2,n) = G2inf(1,1,n);
    G2inf(3,3,n) = (1/2).*a.^(-3).*c.^(-1).*kd(n).^(-2).*(sqrt(-1).*exp(1).^( ...
  sqrt(-1).*a.*c.*kd(n)).*(sqrt(-1)+a.*c.*kd(n))) * coe;
    Gtot(:,:,n) = G0l(:,:,n) + G2inf(:,:,n);
end

%%  infinite integral of coulomb green function
Gcol = zeros(3,3,length(lda));
Gcolxx = zeros(length(lda));
for n = 1:length(lda)
    Gcolxx(n) = 1/4*a.^(-3).*kd(n).^(-2).*(exp(1i.*a.*kd(n))-1i.*a.*kd(n).*(ExpintEi(1i.*a.*kd(n))-1i*pi));
end
for n = 1:length(lda)
Gcol(:,:,n) = [Gcolxx(n),0,0;0,Gcolxx(n),0;0,0,-2*Gcolxx(n)]*coe;
end

%% polarizabiltiy for mat. 
alphaeff = zeros(3,3,length(lda));
alphaE = zeros(3,3,length(lda));
for kk = 1:length(lda)
alphaeff(:,:,kk) = vpa([epsd(kk) - 1,0,0;0,epsd(kk)-1,0;0,0,2.9^2-1]);
end
%% calculate area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%square
if coe==1
[x, y] = meshgrid(0:(Narray-1), 0:(Narray-1));
pos = zeros(2, Narray, Narray); % 2 for (x, y)
V = a^2 * L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%hexagonal
elseif coe==2/sqrt(3)
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
if quarter == 1
valid_points = (element(:,1) >= -1e-5) & (element(:,2) >= -1e-5);
x = element(valid_points, 1);
y = element(valid_points, 2);
V = 0.25 * spacing^2 *sqrt(3)/2 * L;
else
x = element(:,1);
y = element(:,2);
end
V = spacing^2 *sqrt(3)/2 * L;
end

element1 = spacing*[x(:)+shift/spacing,y(:)];

N = size(element1,1);
Ntot(xx) = N;
Scal(xx) = ((floor(Narray*2/sqrt(3))-1)*spacing)^2 *  (3 * sqrt(3) / 2);
corRec = generate_rounded_rectangle(10*1e-9,120*1e-9,40*1e-9);

figure(1)
scatter(element1(:,1),element1(:,2), 5, 'b', 'filled');
hold on
plot(corRec(1,:),corRec(2,:),'r-')
xlim([0, 40e-8])
ylim([0, 40e-8])
% axis equal

%% polarizability for dipole
Pinf = zeros(3,length(lda)); % dipole moment for infinte identical dipoles array
sigma = zeros(1, length(lda)); % absorption for single Pinf (unit W)
for kk = 1:length(lda)
alphaE(:,:,kk) = vpa(eye(3)/ (1 / V * eye(3)/alphaeff(:,:,kk) + k0(kk)^2 * Gcol(:,:,kk)));
end
for kk = 1:length(lda)
alphaE1(:,:,kk) =   Amplifier * vpa(eye(3)/ (1 / V * eye(3)/alphaeff(:,:,kk) + k0(kk)^2 * Gcol(:,:,kk)));
end
 % Qabs = 0.5w Im(p* E)
 % Qext = 0.5w Im(p E0*)
%% polarizability for Si Phy. Rev. B 82, 045404 2010
% n_silver = load('Si n.txt');
% k_silver = load('Si k.txt');
% n_interp = interp1(n_silver(:,1)*1e-6, n_silver(:,2), lda, 'linear','extrap');
% k_interp = interp1(k_silver(:,1)*1e-6, k_silver(:,2), lda, 'linear','extrap');
% epsR = n_interp.^2 - k_interp.^2;
% epsI = 2 * n_interp .* k_interp;
% epsSi = epsR + 1i*epsI;
% nSi = n_interp + 1i*k_interp;
% R = 65*1e-9;
% alphaE = zeros(3,3,length(lda));
% for ii = 1:length(lda)
% alphaE(:,:,ii) = 1i * 6*pi/(k0(ii))^3 * eps0 * calculate_a1(nSi(ii), k0(ii)*R) *eye(3);
% end
%% polarizability for Silver J.Phys.:Photonics1(2019)015004
% n_silver = load('JohnSonSilver n.txt');
% k_silver = load('JohnSonSilver k.txt');
% n_interp = interp1(n_silver(:,1)*1e-6, n_silver(:,2), lda, 'linear','extrap');
% k_interp = interp1(k_silver(:,1)*1e-6, k_silver(:,2), lda, 'linear','extrap');
% epsR = n_interp.^2 - k_interp.^2;
% epsI = 2 * n_interp .* k_interp;
% epsAg = epsR + 1i*epsI;
% nAg = n_interp + 1i*k_interp;
% R = 50*1e-9;
% 
% for ii = 1:length(lda)
% alphaE(:,:,ii) = 3/2 / (k0(ii))^3 * eye(3) * calculate_t1E(R, lda(ii), 1, nAg(ii));
% end
%% Plasmonic Green function
% en function
omegac = 1.7712; % eV resonance of plasma 1.734 1.7712
gammac = 0.2031; % eV decay rate of plasma 0.1946 0.2031/2
pos = [element1, 0.5*1e-9*ones(N, 1)] * 1e9 + shift;

maxE = mphmax(model,'ewfd2.normE',3,'selection',23);
                    Ex = mphinterp(model,'ewfd2.Ex-ewfd2.Ebx','coord', pos');
                    Ey = mphinterp(model,'ewfd2.Ey-ewfd2.Eby','coord', pos');
                    Ez = mphinterp(model,'ewfd2.Ez-ewfd2.Ebz','coord', pos');
E0 = [Ex;Ey;Ez];
u = E0/maxE;

tic;

lorentzShape = 1./(omegac^2-eV.^2-1i*gammac.*eV);
[maxValue, maxIndex] = max(lorentzShape);
Ptot = zeros(3*N,length(lda));
Pd = zeros(3*N,length(lda));
flag = false;
end
%% Phase of dipoles
Phase = zeros(N,N);
upperTri = 2*pi*(rand(N) + 1i * rand(N));  % chi(r1,r2) = conj(chi(r2,r1));                          
upperTri = triu(upperTri, 1); 
phase = upperTri + upperTri';
phase = phase + eye(N);
if strcmp(current_config.param_name , 'Repeating')
phase_tensor = reshape(exp(phase),[1,1,N,N]);
phase_tensor = bsxfun(@times,phase_tensor,reshape(eye(3),[3,3,1,1]));
else
phase_tensor = reshape(exp(zeros(N,N)),[1,1,N,N]);
phase_tensor = bsxfun(@times,phase_tensor,reshape(eye(3),[3,3,1,1]));
end
Egnr = zeros(3,N,length(lda));
%% Core calculation area
parfor kk = 1:length(lda)
M = ones(3,3,N,N);
Mdip = ones(3,3,N,N);
Gmnp = zeros(3,3,N,N);
Gij = zeros(3,3,N,N);
Gfinal = zeros(3,3,N,N);
[iIndex, jIndex] = ndgrid(1:N, 1:N);
maskOffDiag = (iIndex ~= jIndex);  
maskDiag    = (iIndex == jIndex);

    [i, j] = ndgrid(1:N, 1:N); 
    p1 = element1(i(:), :); 
    p2 = element1(j(:), :); 
    x1 = p1(:, 1); 
    y1 = p1(:, 2); 
    z1 = zeros(size(x1)); 
    x2 = p2(:, 1); 
    y2 = p2(:, 2);
    z2 = zeros(size(x2));
    Gij = tensorG_vector(kd(kk), x1, y1, z1, x2, y2, z2); % dipole-dipole Green function
    Gij = pagemtimes(reshape(Gij, 3, 3, N, N),phase_tensor); 
    Gmnp = tensorG_mnp(u(:,i),u(:,j),N) * ( A/(omegac^2-eV(kk)^2-1i*gammac*eV(kk)) ); % plasmonic Green function
    Gmnp = Amplifier * pagemtimes(reshape(Gmnp,3,3,N,N),phase_tensor);
    Gfinal = Gij + Gmnp;
    alphaEmatrix = repmat(alphaE(:,:,kk),1,1,N,N);
    M = -k0(kk)^2 * (pagemtimes(alphaEmatrix, Gij) + pagemtimes(alphaEmatrix, Gmnp)) ; 
    Mdip = -k0(kk)^2 * (pagemtimes(alphaEmatrix, Gij));
 
diagMask = false(N, N); 
diagMask(1:N+1:end) = true; 
diagMask = repmat(reshape(diagMask, [1, 1, N, N]), [3, 3, 1, 1]);

% Mdip(:,:,diagIndices,diagIndices) = eye(3)
eyeBlock = repmat(eye(3), [1, 1, N]); 
Mdip(diagMask) = eyeBlock(:); 

% M(:,:,diagIndices,diagIndices) = eye(3) + M(:,:,diagIndices,diagIndices)
M(diagMask) = M(diagMask) + eyeBlock(:);

M_flat = reshape(permute(M, [1, 3, 2, 4]), [3*N, 3*N]); 
Mdip_flat = reshape(permute(Mdip, [1, 3, 2, 4]), [3*N, 3*N]);

W = M_flat;
Wdip = Mdip_flat;

Ufreeiso = eps0 * alphaE(:,:,kk) * E0free;
% Uiso = alphaE(kk)*eye(3) * E0;
Ufree = repmat(Ufreeiso, N, 1);
U = zeros(3*N,1);

E0_freq = E0 * ( gammac*omegac/(omegac^2-eV(kk)^2-1i*gammac*eV(kk)) );
% E0_freq = E0;
U_segments = pagemtimes(eps0 * alphaE(:,:,kk), reshape(E0_freq,3,1,N));
U = reshape(U_segments,[],1);

Ptot(:,kk) = eye(3*N) *  pinv(W) * U;
Pd(:,kk) = eye(3*N)/(Wdip)*Ufree;
Pd_reshaped = reshape(Pd(:,kk), [3, N]);
Ptot_reshaped = reshape(Ptot(:,kk), [3, N]);

Egnr_kk = zeros(3,N);
Edipgnr_kk = zeros(3,N);
    for ii = 1:N
        for jj = 1:N
        Egnr_kk(:,ii) = Egnr_kk(:,ii) + k0(kk)^2 / eps0 * Gfinal(:,:,jj,ii) * Ptot_reshaped(:,jj);
        Edipgnr_kk(:,ii) = Edipgnr_kk(:,ii) + k0(kk)^2 / eps0 * Gij(:,:,jj,ii) * Pd_reshaped(:,jj);
        end
    end
if quarter == 1
    frac = 4;
else
    frac = 1;
end
Qext_tot(xx,kk) = frac * omega0(kk) * imag(sum(sum(conj(E0 * ( gammac * omegac/(omegac^2-eV(kk)^2-1i*gammac*eV(kk)) )) .* (Ptot_reshaped), 1)));
Qext_dip(xx,kk) = frac * omega0(kk) * imag(sum(sum(conj(E0free) .* (Pd_reshaped), 1)));
% Qext_dip(xx,kk) = 4*pi*k0(kk) * imag(sum(sum(conj(E0free) .* (Pd_reshaped), 1)));
% Qext_comp(xx,kk) = Qext_comp(xx,kk) + 0.5 * omega0(kk) * real(-1i * sum(sum(conj(E0free) .* (Ptot((ii-1)*3 + (1:3),kk)))));
Qabs_tot(xx,kk) = frac * omega0(kk) * (imag(sum(sum(Ptot_reshaped .*(conj(eye(3)/alphaE(:,:,kk)/eps0)*conj(Ptot_reshaped)), 1)))-sum(sum(2/3*(k0(kk))^3*(abs(Ptot_reshaped)).^2)));
Qabs_dip(xx,kk) = frac * omega0(kk) * (imag(sum(sum(Pd_reshaped .*(conj(eye(3)/alphaE(:,:,kk)/eps0)*conj(Pd_reshaped)), 1)))-sum(sum(2/3*(k0(kk))^3*(abs(Pd_reshaped)).^2)));
% Qabs_tot(xx,kk) = -omega0(kk) * (imag(sum(sum(conj(Ptot_reshaped) .*((eye(3)/alphaE(:,:,kk)/eps0)*conj(Ptot_reshaped)), 1))));
% Qabs_dip(xx,kk) = -omega0(kk) * (imag(sum(sum(conj(Pd_reshaped) .*((eye(3)/alphaE(:,:,kk)/eps0)*conj(Pd_reshaped)), 1))));

%% infinite dipoles
Pinf(:,kk) = eps0 * eye(3)/(eye(3)/alphaE(:,:,kk) - k0(kk)^2 * Gtot(:,:,kk)) * E0free;
sigma(kk) = omega0(kk) * sum(imag((Pinf(:,kk)).*conj(eye(3)/(alphaE(:,:,kk))/eps0 * Pinf(:,kk)))-sum(sum(2/3*(k0(kk))^3*(abs(Pinf(:,kk))).^2)));
end
timeElapsed = toc;
fprintf('calculation time: %.4f sec\n', timeElapsed);
Qabs_tot_avg(xx,:) = sum(Qabs_tot,1)/xx;
Qext_tot_avg(xx,:) = sum(Qext_tot,1)/xx;
Qabs_dip_avg(xx,:) = sum(Qabs_dip,1)/xx;
end
%% Save
% Close the waitbar after completion
results.lda = lda;
results.Qavg = Qext_tot_avg;
results.Qdipavg = Qabs_dip_avg;  
results.Qabsavg = Qabs_tot_avg;
% results.Qscaavg = Qext_dip_avg(size(Qext_dip_avg,1),:)-Qabs_tot_avg(size(Qabs_tot_avg,1)-1);
Qsca = Qext_tot - Qabs_tot;
Qsca_avg = Qext_tot_avg - Qabs_tot_avg;
results.Qext = Qext_tot;
results.Qdip = Qext_dip;
results.Qabs = Qabs_tot;
results.Qsca = Qext_tot - Qabs_tot;
results.Qdipabs = Qabs_dip;
results.Alpha = alphaE;
results.alpha2D = alphaeff;
results.ImGmnp = ImGmnpxx;
results.ImGij = ImGijxx;
results.ReGmnp = ReGmnpxx;
results.ReGij = ReGijxx;
results.Pi = Ptot;
results.Qinf = sigma;
results.Qcomp = Qext_comp;
results.Ntot = Ntot;
results.para = paramArray;
FIGS = strcat('WSe2na1.15_145Kspacing',num2str(spacing*1e9),'nmN=',num2str((Narray)),'Amax=',num2str(max(paramArray)));
% filename = strcat('WSe2na1.15_145Kspacing',num2str(spacing*1e9),'nmN=',num2str((Narray)),'Repeatmax=',num2str(max(paramArray)),'.mat');
% save(filename,'results');
clear all
clc;
% Choose ur parameter

%% constant
eps0 = 8.8541878188*1e-12; %8.8541878188*1e-12 F⋅m−1 vacuum dielectric constant
econst = 1.602176634 * 1e-19; %C
hbar = 1.054571817*1e-34 ; %J`s
cconst = 299792458;
%% parameters
% dipoles properties
radius = 9 * 1e-9;
initialphase = 0; %initial phase of dipoles
theta = 0 * pi/180;
alpha = 1; %coefficient for inhomogeneous frequency distribution
L = 1 * 1e-9; % height of 2D mat.
Narray = 15; % element number generator
spacing = 5*1e-9; %distance between dipoles
fA1 = 0.9; %1.9[eV]
omegaA1 = 1.6674; %2.014[eV]  1.6674 1.7072
gammaA1 = 0.035; % 0.0281 0.018
epsinf = 1; 
n1 = 1; % refractive index of air
Amplifier = 1; % amplitude for polarizability
shiftx = 0;
shifty = 0;
%postion of the origin point of array
A = 1.79e6; % coefficient for Gmnp(R,R0) 1.79e6
flag = true;
configs = struct(...
    'use_R', {1, 2, 3, 4, 5, 6, 7, 8}, ...
    'param_name', {'spacing', 'Narray', 'radius','Repeating','sigma_Gaus','alpha','theta','initialphase'}, ...
    'paramArray', {linspace(1,5,20)*1e-9, ...
    [5,8,10], ...
    linspace(10,100,10), ...
    [1:1:5], ...
    linspace(10,100,10) ...
    linspace(0,1,11), ...
    0*pi/180 ...
    linspace(0,pi,5)});
% alpha is for the combination between w = a w+ + (1-a) w-
use_R = 1; 
lda = linspace(700,800,200) * 1e-9; % incident wavelength
% lda = linspace(400,700,100) * 1e-9;
E0free = [sqrt(1.3707e17);0;0];%V/m amplitude of incident electric field in free space
sigma_Gaus = 1e10; %Standard deviation of Gaussian distribution

coe = 2/sqrt(3); % 1 for square, 2/sqrt(3) for Hexago
quarter = 0;

k0 = 2*pi./lda;
kd = n1*k0;
omega0 = 299792458 * k0;
eV = 299792458 * 6.58211899e-16 * (2 * pi) ./ lda;

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
eval([current_config.param_name '=' num2str(paramArray(xx))]);
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
%% Just for fun
if strcmp(current_config.param_name , 'theta')

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

% 
% 提取 x 和 y 坐标
% x = element1(:, 1);
% y = element1(:, 2);

% 定义旋转矩阵
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% 对每个点进行旋转
rotated_points = (R * [x'; y'])';

% 更新 element1 为旋转后的坐标
element0 = rotated_points;
element0full = [element0, L*ones(size(element0,1),1)];
% element0full = [element0, zeros(size(element0,1),1)];
element1full = [element1, zeros(size(element0,1),1)];
element1 = unique([element0full; element1full], 'rows');
Scal(xx) = ((floor(Narray*2/sqrt(3))-1)*spacing)^2 *  (3 * sqrt(3) / 2);
corRec = generate_rounded_rectangle(10*1e-9,120*1e-9,40*1e-9);
x = element1(:,1);
y = element1(:,2);


%% calculate area
else
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
x = element(:, 1);
y = element(:, 2);
end
V = spacing^2 *sqrt(3)/2 * L;
end

% element1 = spacing*[x(:)+shiftx/spacing, y(:)+shifty/spacing];
% center = [shiftx, shifty];
% % 计算相对于中心的坐标
% dx = element1(:,1) - center(1);
% dy = element1(:,2) - center(2);
% abs_dx = abs(dx);
% abs_dy = abs(dy);
% % 计算正六边形允许的最大y值
% max_dy_allowed = sqrt(3) * min(radius/2, radius - abs_dx);
% % 判断点是否在正六边形内
% inside_hexagon = (abs_dx <= radius) & (abs_dy <= max_dy_allowed);
% element_inside = [element1(inside_hexagon, 1), element1(inside_hexagon, 2)];
% element1 = [element_inside, zeros(size(element_inside, 1), 1)];

element1 = spacing*[x(:)+shiftx/spacing,y(:)+shifty/spacing];
center = [shiftx,shifty];
dist_from_center = sqrt((element1(:,1) - center(1)).^2 + (element1(:,2) - center(2)).^2);
inside_circle = dist_from_center <= radius;
element_inside = [element1(inside_circle,1),element1(inside_circle,2)];
element1 = [element_inside,zeros(size(element_inside,1),1)];

% element1 = [element1, zeros(size(element1,1),1)];
Scal(xx) = ((floor(Narray*2/sqrt(3))-1)*spacing)^2 *  (3 * sqrt(3) / 2);
corRec = generate_rounded_rectangle(10*1e-9,120*1e-9,40*1e-9);
end

figure
scatter(element1(:,1),element1(:,2), 5, 'b', 'filled');
tta = linspace(0, 2*pi, 100);  % 生成角度数据
x_circle = radius * cos(tta);  % 圆框的 x 坐标
y_circle = radius * sin(tta);  % 圆框的 y 坐标
hold on
plot(x_circle, y_circle, 'r', 'LineWidth', 2); 
xlim([-radius*1.2,radius*1.2]);
ylim([-radius*1.2,radius*1.2]);
axis equal
outputFolder = strcat('E:\ljx\记录\coupledipolemodel\radius\radius_', num2str(radius), 'nm');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder); 
end
    fileName = fullfile(outputFolder, ['gap_',num2str(spacing),'.jpg']);  
    saveas(gcf, fileName, 'jpg');  

%% Area disorder

% radius = 10*1e-9;
% 
% random_index = randi(size(element1, 1)); 
% origin = element1(random_index, :); 
% 
% distances = sqrt((element1(:,1) - origin(1)).^2 + (element1(:,2) - origin(2)).^2);
% 
% inside_region = distances <= radius;
% 
% element1(inside_region, :) = [];
%
% 
% figure
% scatter(element1(:,1),element1(:,2), 5, 'b', 'filled');
% axis equal
% hold on

N = size(element1,1);
Ntot(xx) = N;
%% polarizabiltiy for mat. 
r_sq = x.^2 + y.^2;
omegaA2 = omegaA1 * (exp(-r_sq/(2*sigma_Gaus^2)*0));
omegaA3 = omegaA1 + omegaA1 * (1 - exp(-r_sq/(2*sigma_Gaus^2)));
omegaA0 = alpha * omegaA2 + (1-alpha) * omegaA3;

% figure(1)
% scatter3(x(:), y(:), omegaA0(:), 50, omegaA0(:), 'filled') 
% colormap(jet)
% % colorbar20
% title('Resonance Frequency Distribution')
% xlabel('x'), ylabel('y'), zlabel('\omega_0^\prime')
% view(-30, 45) 
% grid on
% z_max = omegaA1;
% z_min = 1.53;
% zlim([z_min, z_max])

epsd = epsinf + fA1./(omegaA0.^2-eV.^2-1i*gammaA1.*eV); % dielectric function for dipole
diag_components = cat(3, epsd - 1, epsd - 1, repmat(2.9^2 - 1, size(epsd)));
diag_components = permute(diag_components, [3 1 2]);
alphaeff = zeros(3, 3, N, length(lda));

component_index = [1 2 3]; % 对应xx,yy,zz分量

for k = 1:length(lda)
    for n = 1:Ntot(xx)
        diag_values = diag_components(:, n, k); % 提取当前格点的三个分量 [3x1]
        alphaeff(:,:,n,k) = diag(diag_values);  % 生成3x3对角矩阵
    end
end
deltaMax(xx) = max(omegaA0) - min(omegaA0);
%% polarizability for dipole
Pinf = zeros(3,length(lda)); % dipole moment for infinte identical dipoles array
sigma = zeros(1, length(lda)); % absorption for single Pinf (unit W)

alphaE = zeros(3, 3, N, length(lda)); % 预分配 alphaE (3x3xNxn)

% 主循环
for kk = 1:length(lda)
    % 提取当前波长的有效极化率
    current_alpha = alphaeff(:,:,:,kk); % 3x3xN
    
    % 计算 term = (current_alpha)^{-1} / V + k0^2 * Gcol
    term = pagemtimes(pagepinv(current_alpha), 1/V) + ...
           k0(kk)^2 * Gcol(:,:,kk); % 假设 Gcol 是 3x3xNxn
    
    % 计算 alphaE = (eye(3) + term)^{-1}
    alphaE(:,:,:,kk) = pagepinv(term);
end

 % Qabs = 0.5w Im(p* E)
 % Qext = 0.5w Im(p E0*)
tic;

Pd = zeros(3*N,length(lda));
flag = false;
end

%% Phase of dipoles
Phase = zeros(N,N);
upperTri = 0 * 2*pi*(rand(N) + 1i * rand(N));  % chi(r1,r2) = conj(chi(r2,r1));                          
upperTri = triu(upperTri, 1); 
phase = upperTri + upperTri';
phase = phase + eye(N);
if strcmp(current_config.param_name , 'Repeating')
phase_tensor = reshape(exp(phase),[1,1,N,N]);
phase_tensor = bsxfun(@times,phase_tensor,reshape(eye(3),[3,3,1,1]));
else
    phase_tensor = reshape(exp(1i*initialphase*ones(N,N)),[1,1,N,N]);
    phase_tensor = bsxfun(@times,phase_tensor,reshape(eye(3),[3,3,1,1]));
end
for i = 1:N
    phase_tensor(:,:,i,i) = eye(3);
end
Egnr = zeros(3,N,length(lda));
%% Core calculation area
parfor kk = 1:length(lda)
Mdip = ones(3,3,N,N);
Gij = zeros(3,3,N,N);
[iIndex, jIndex] = ndgrid(1:N, 1:N);
maskOffDiag = (iIndex ~= jIndex);  
maskDiag    = (iIndex == jIndex);

    [i, j] = ndgrid(1:N, 1:N); 
    p1 = element1(i(:), :); 
    p2 = element1(j(:), :); 
    x1 = p1(:, 1); 
    y1 = p1(:, 2); 
    z1 = p1(:, 3); 
    x2 = p2(:, 1); 
    y2 = p2(:, 2);
    z2 = p2(:, 2);
    Gij = tensorG_vector(kd(kk), x1, y1, z1, x2, y2, z2); % dipole-dipole Green function
    Gij = pagemtimes(reshape(Gij, 3, 3, N, N),phase_tensor); 

    Mdip = -k0(kk)^2 * (pagemtimes(alphaE(:,:,:,kk), Gij));
 
diagMask = false(N, N); 
diagMask(1:N+1:end) = true; 
diagMask = repmat(reshape(diagMask, [1, 1, N, N]), [3, 3, 1, 1]);

% Mdip(:,:,diagIndices,diagIndices) = eye(3)
eyeBlock = repmat(eye(3), [1, 1, N]); 
Mdip(diagMask) = eyeBlock(:); 

Mdip_flat = reshape(permute(Mdip, [1, 3, 2, 4]), [3*N, 3*N]);

Wdip = Mdip_flat;
E0freeRepeat = repmat(E0free,1,Ntot(xx));
E0freeRepeat = reshape(E0freeRepeat, [3, 1, N]);
Ufree = squeeze(eps0 * pagemtimes(alphaE(:,:,:,kk) , E0freeRepeat));
Ufree_reshaped = reshape(Ufree, [3*N, 1]);
Pd(:,kk) = eye(3*N)/(Wdip)*Ufree_reshaped;
Pd_reshaped = reshape(Pd(:,kk), [3, N]);

Edipgnr_kk = zeros(3,N);
    for ii = 1:N
        for jj = 1:N
        Edipgnr_kk(:,ii) = Edipgnr_kk(:,ii) + k0(kk)^2 / eps0 * Gij(:,:,jj,ii) * Pd_reshaped(:,jj);
        end
    end
if quarter == 1
    frac = 4;
else
    frac = 1;
end
Qext_dip(xx,kk) = frac * omega0(kk) * imag(sum(sum(conj(E0free) .* (Pd_reshaped), 1)));
Qabs_dip(xx,kk) = frac * omega0(kk) * (imag(sum(sum(Pd_reshaped .*squeeze(pagemtimes(conj(pagepinv(alphaE(:,:,:,kk))/eps0) , conj(reshape(Pd_reshaped,3,1,N)))), 1)))-sum(sum(2/3*(k0(kk))^3*(abs(Pd_reshaped)).^2)));

%% infinite dipoles
% Pinf(:,kk) = eps0 * eye(3)/(eye(3)/alphaE0(:,:,kk) - k0(kk)^2 * Gtot(:,:,kk)) * E0free;
% sigma(kk) = omega0(kk) * sum(imag((Pinf(:,kk)).*conj(eye(3)/(alphaE0(:,:,kk))/eps0 * Pinf(:,kk)))-sum(sum(2/3*(k0(kk))^3*(abs(Pinf(:,kk))).^2)));
end
timeElapsed = toc;

fprintf('calculation time: %.4f sec\n', timeElapsed);
Qabs_dip_avg(xx,:) = sum(Qabs_dip,1)/xx;
end
%% Save
% Close the waitbar after completion
results.lda = lda;
results.Qabs = Qabs_dip;
results.Qdipavg = Qabs_dip_avg;  
results.Qdip = Qext_dip;
results.Alpha = alphaE;
results.alpha2D = alphaeff;
results.ImGmnp = ImGmnpxx;
results.ImGij = ImGijxx;
results.ReGmnp = ReGmnpxx;
results.ReGij = ReGijxx;
% results.Qinf = sigma;
results.Qcomp = Qext_comp;
results.Ntot = Ntot;
results.para = paramArray;
FIGS = strcat('WSe2na1.15_145Kspacing',num2str(spacing*1e9),'nmN=',num2str((Narray)),'Amax=',num2str(max(paramArray)));
% filename = strcat('WSe2na1.15_145Kspacing',num2str(spacing*1e9),'nmN=',num2str((Narray)),'Repeatmax=',num2str(max(paramArray)),'.mat');
% save(filename,'results');
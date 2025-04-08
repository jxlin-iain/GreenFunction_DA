% This code contains Green function convergence test (inculding basic integral area; different part of green function; limitation of green function) You need to take use of it by some function code.
% The calculation time may be long (takes hours) when the NarryPara increased to 3e2 and
% more. The parameters for the final part (different starting point) should
% be re-define before running all.
clear all
clc
lda = linspace(600,900,20)*1e-9; % incident wavelength
% lda = 700 * 1e-9;
k0 = 2 * pi ./ lda;  % wave vevtor in free space
eV = 299792458 * 6.58211899e-16 * (2 * pi) ./ lda; %[eV]
spacing = 5*1e-9; % distance between elements
% NarrayPara = [10,50,60,70,80,90,1e2,2e2,3e2,4e2,5e2]; % element number generator
NarrayPara = [10,50,60,70,80,90,1e2]; % element number generator
% NarrayPara = [4];
% startPoint = [10,20,30,40,50,60,70];
startPoint = 1;
c = startPoint; % lower limitation of integral / default: 1
a = spacing;
b = NarrayPara; 
G0 = zeros(3,3,length(NarrayPara),length(lda));
G0far = zeros(3,3,length(NarrayPara),length(lda));
G0col = zeros(3,3,length(NarrayPara),length(lda));
G0colRe = zeros(3,3,length(NarrayPara),length(lda));
coe = 2/sqrt(3); % the ratio between int & sum, i.e. sum=int*coe

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare between integral of 1 & calculation area display
int = (-c+b.^2)*pi; % integral of rdrdthtea
for mm = 1:length(NarrayPara)
Narray = NarrayPara(mm);
%% Hexagonal array generation

element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);
 mask = R < Narray;
x = x(mask);
y = y(mask);
element = a*[x(:),y(:)];
ATransposed = element';

[uniqueRows, ~, idx] = unique(ATransposed, 'rows', 'stable');
if mm==1
if size(uniqueRows, 1) < size(ATransposed, 1)
    disp('There are duplicate 2x1 vectors in the matrix.');
else
    disp('No duplicate position found.');
end
end
origin_idx = find(element(:,1) == 0 & element(:,2) == 0, 1);

if ~isempty(origin_idx)
    element = circshift(element, -(origin_idx - 1));
end
% scatter(element(:,1),element(:,2))
element_number(mm) = size(element, 1);  % Total iterations for the inner loop

end
figure(1)
scatter(element(:,1),element(:,2))
title('reduced element structure (please zoom in)')
figure(2)
plot(b,int*coe)
hold on
plot(b,element_number,'--'); % coefficient sqrt(3)/2 has to be taken into account for the geometry reason
ylabel('unit summation/element number')
xlabel('normalized radius')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');
title('comparision of area between integral and summation')
saveas(gcf, 'unit.png');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare between integral of G~1/r
Gfar = zeros(length(b),length(lda));
for n = 1:length(lda)
Gfar(:,n) = (sqrt(-1)*(1/4)).*a.^(-2).*(exp(1).^(sqrt(-1).*a.*k0(n))+(-1).*exp(1) ...
  .^(sqrt(-1).*a.*b.*k0(n))).*k0(n).^(-1);
end
for mm = 1:length(NarrayPara)
Narray = NarrayPara(mm);
%% Hexagonal array generation
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);

% mask = (R < Narray) & (R > startPoint(mm));
 mask = R < Narray;
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
        Gij = tensorGfar(k0(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0far(:,:,mm,n) = G0far(:,:,mm,n) + Gij;
    end
end
end
figure(3)
plot(element_number,real(Gfar(:,1))*coe)
hold on
plot(element_number,squeeze(real(G0far(1,1,:,1))),'--')
xlabel('element number')
ylabel('sumGfar')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  title('real part of Gfar between integral and summation')
saveas(gcf, 'ReGfar.png');
figure(4)
plot(element_number,imag(Gfar(:,1))*coe)
hold on
plot(element_number,squeeze(imag(G0far(1,1,:,1))),'--')
xlabel('element number')
ylabel('sumGfar')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none'); 
title('Imaginary part of Gfar between integral and summation')
saveas(gcf, 'ImGfar.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare between integral of coulomb Green function
% Gcol has already integral by theta direction
rmin = 1;
for ii = 1:length(lda)
for n = 1:length(b)
r = rmin:1:b(n);
Gcol(n,ii) = sum((1/4).*a.^(-3).*exp(1).^(sqrt(-1).*a.*k0(ii).*r).*k0(ii).^(-2).*r.^(-2)); % approximate the integral by discrete summation
end
end
for mm = 1:length(NarrayPara)
Narray = NarrayPara(mm);
%% Hexagonal array generation
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);

% mask = (R < Narray) & (R > startPoint(mm));
 mask = R < Narray;
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
        Gij = tensorG_coulomb(k0(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0col(:,:,mm,n) = G0col(:,:,mm,n) + Gij;
    end
end
end
%1.0672
figure(5)
plot(element_number,real(Gcol(:,1))*coe)
hold on
plot(element_number,squeeze(real(G0col(1,1,:,1))),'--')
xlabel('element number')
ylabel('sumGcol')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  
title('real part of Gcol between integral and summation')
saveas(gcf, 'ReGcol.png');
figure(6)
plot(element_number,imag(Gcol(:,1))*coe)
hold on
plot(element_number,squeeze(imag(G0col(1,1,:,1))),'--')
xlabel('element number')
ylabel('sumGcol')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  
title('Imaginary part of Gcol between integral and summation')
saveas(gcf, 'ImGcol.png');
% real(squeeze(G0col(1,1,:,1)))./real(Gcol(:,1));
% imag(squeeze(G0col(1,1,:,1)))./real(Gcol(:,1));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare between integral of coulomb Green function without exp

for n = 1:length(lda)
GcolRe(:,n) = (1/4).*a.^(-3).*((-1)+b).*b.^(-1).*k0(n).^(-2);
end
for mm = 1:length(NarrayPara)
Narray = NarrayPara(mm);
%% Hexagonal array generation
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);

% mask = (R < Narray) & (R > startPoint(mm));
 mask = R < Narray;
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
        Gij = tensorG_coulombRe(k0(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0colRe(:,:,mm,n) = G0colRe(:,:,mm,n) + Gij;
    end
end
end

figure(7)
plot(element_number,real(GcolRe(:,1))*coe)
hold on
plot(element_number,squeeze(real(G0colRe(1,1,:,1))),'--')
xlabel('element number')
ylabel('sumGcolRe')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  title('real part of GcolRe between integral and summation')
saveas(gcf, 'GcolRe.png');

% real(squeeze(G0colRe(1,1,:,1)))./real(GcolRe(:,1))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare between integral of full Green function
for n = 1:length(lda)
G2x(:,n) = ((1/4).*a.^(-3).*b.^(-1).*c.^(-1).*k0(n).^(-2).*((sqrt(-1)*(-1)).*c.* ...
  exp(1).^(sqrt(-1).*a.*b.*k0(n)).*((sqrt(-1)*(-1))+a.*b.*k0(n))+b.*exp(1) ...
  .^(sqrt(-1).*a.*c.*k0(n)).*(1+sqrt(-1).*a.*c.*k0(n))));
G2z(:,n) = (1/2).*a.^(-3).*c.^(-1).*k0(n).^(-2).*(sqrt(-1).*exp(1).^( ...
  sqrt(-1).*a.*c.*k0(n)).*(sqrt(-1)+a.*c.*k0(n)));
  % % infinete integral
end
for mm = 1:length(NarrayPara)
Narray = NarrayPara(mm);
%% Hexagonal array generation
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);

% mask = (R < Narray) & (R > startPoint(mm));
 mask = R < Narray;
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
        Gij = tensorG(k0(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0(:,:,mm,n) = G0(:,:,mm,n) + Gij;
    end
end
end

figure(8)
plot(element_number,real(G2x(:,1))*(coe))
hold on
plot(element_number,squeeze(real(G0(1,1,:,1))),'--')
xlabel('element number')
ylabel('sumGcol')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  
title('real part of Gfullxx between integral and summation')
saveas(gcf, 'ReGfullxx.png');
figure(9)
plot(element_number,imag(G2x(:,1))*(coe))
hold on
plot(element_number,squeeze(imag(G0(1,1,:,1))),'--')
xlabel('element number')
ylabel('sumGcol')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  
title('Imaginary part of Gfullxx between integral and summation')
saveas(gcf, 'ImGfullxx.png');

figure(10)
plot(element_number,real(G2z(:,1))*(coe))
hold on
plot(element_number,squeeze(real(G0(3,3,:,1))),'--')
xlabel('element number')
ylabel('sumGcol')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  
title('real part of Gfullzz between integral and summation')
saveas(gcf, 'ReGfullzz.png');
figure(11)
plot(element_number,imag(G2z(:,1))*(coe))
hold on
plot(element_number,squeeze(imag(G0(3,3,:,1))),'--')
xlabel('element number')
ylabel('sumGcol')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  
title('Imaginary part of Gfullzz between integral and summation')
saveas(gcf, 'ImGfullzz.png');

% real(squeeze(G0(1,1,:,1)))./real(G2(:,1))
% imag(squeeze(G0(1,1,:,1)))./imag(G2(:,1))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% different starting point

startPoint = 0.01:0.04:0.8; 
% startPoint = 1;
NarrayPara = 1e2;
c = floor(startPoint.*NarrayPara); % lower limitation of integral, corresponding to m in the pdf.
b = NarrayPara;
G0s = zeros(3,3,length(startPoint),length(lda));
G2x=zeros(length(c),length(lda));
for n = 1:length(lda)
G2x(:,n) = ((1/4).*a.^(-3).*b.^(-1).*c.^(-1).*k0(n).^(-2).*((sqrt(-1)*(-1)).*c.* ...
  exp(1).^(sqrt(-1).*a.*b.*k0(n)).*((sqrt(-1)*(-1))+a.*b.*k0(n))+b.*exp(1) ...
  .^(sqrt(-1).*a.*c.*k0(n)).*(1+sqrt(-1).*a.*c.*k0(n))));
% G2lzz(:,n) = (1/2).*a.^(-3).*c.^(-1).*k0(n).^(-2).*(sqrt(-1).*exp(1).^( ...
%   sqrt(-1).*a.*c.*k0(n)).*(sqrt(-1)+a.*c.*k0(n)));
% G2(n) =  ((1/4).*a.^(-3).^(-1).*c.^(-1).*k0(n).^(-2).*(exp(1) ...
  % .^(sqrt(-1).*a.*c.*k0(n)).*(1+sqrt(-1).*a.*c.*k0(n)))) * (2/sqrt(3))^3;
end
for mm = 1:length(startPoint)
Narray = NarrayPara;
%% Hexagonal array generation
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);

mask = (R < Narray) & (R > c(mm));
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
        Gij = tensorG(k0(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0s(:,:,mm,n) = G0s(:,:,mm,n) + Gij;
    end
end
end
figure(10)
scatter(element(:,1),element(:,2))
title('reduced element structure')
saveas(gcf, 'IntArea.png');
figure(11)
plot(c,real(G2x(:,1))*(coe))
hold on
plot(c,squeeze(real(G0s(1,1,:,1))),'--')
ylabel('sumGcol')
xlabel('starting normalized radius')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  
title('real part of Gfull between integral and summation')
saveas(gcf, 'ReGfulls.png');
figure(12)
plot(c,imag(G2x(:,1))*(coe))
hold on
plot(c,squeeze(imag(G0s(1,1,:,1))),'--')
ylabel('sumGcol')
xlabel('starting normalized radius')
hLegend = legend('int','sum');
set(hLegend, 'Box', 'off');       
set(hLegend, 'Color', 'none');  
title('Imaginary part of Gfull between integral and summation')
saveas(gcf, 'ImGfulls.png');
 % save('differenstartingpoint.mat', 'G0s');

 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% summation + integral vs summation (cacluated area comparison)
stopPoint = 0.01:0.04:0.8; 
NarrayPara = 1e2;
%%%% right hand side
element = generate_hexagonal_structure(floor(NarrayPara*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);
mask = R < NarrayPara;
x = x(mask);
y = y(mask);
element = a*[x(:),y(:)];
origin_idx = find(element(:,1) == 0 & element(:,2) == 0, 1);

if ~isempty(origin_idx)
    element = circshift(element, -(origin_idx - 1));
end
element_numberR = size(element, 1);  % Total iterations for the inner loop

%%%% left hand-side
c = stopPoint*NarrayPara;
b = NarrayPara;
element_numberL = zeros(1,length(c));
intL = (-c.^2+b.^2)*pi; % integral of rdrdthtea
for mm = 1:length(stopPoint)
    Narray = stopPoint(mm)*NarrayPara;
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);
mask = R < Narray;
 % mask = R < Narray;
x = x(mask);
y = y(mask);
element = a*[x(:),y(:)];
origin_idx = find(element(:,1) == 0 & element(:,2) == 0, 1);

if ~isempty(origin_idx)
    element = circshift(element, -(origin_idx - 1));
end
element_numberL(mm) = size(element, 1);  % Total iterations for the inner loop
end
LEFTxx = element_numberL + intL*coe;
RIGHTxx = element_numberR;

diff = RIGHTxx - LEFTxx;

figure(13)
plot(c,LEFTxx)
ylabel('calculatied area/element number')
xlabel('normalized stop summation/start integral radius')
title('sum + int');
saveas(gcf, 'unit_comp.png');

figure(14)
plot(c,element_numberL)
hold on
plot(c,intL,'--')
legend('sum','int')
xlabel('normalized stop summation/start integral radius')
title('sum vs int');
saveas(gcf, 'unit_sumvsint.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% summation + integral vs summation
stopPoint = 0.01:0.04:0.9; 
NarrayPara = 3e2;
G0r = zeros(3,3,length(lda));

%%%% right hand side
element = generate_hexagonal_structure(floor(NarrayPara*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);
mask = R < NarrayPara;
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
        Gij = tensorG(k0(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0r(:,:,n) = G0r(:,:,n) + Gij;
    end
end
%%%% left hand-side
c = stopPoint*NarrayPara;
b = NarrayPara;
G0l = zeros(3,3,length(stopPoint),length(lda));
G2l = zeros(length(c),length(lda));
for n = 1:length(lda)
G2l(:,n) = ((1/4).*a.^(-3).*b.^(-1).*c.^(-1).*k0(n).^(-2).*((sqrt(-1)*(-1)).*c.* ...
  exp(1).^(sqrt(-1).*a.*b.*k0(n)).*((sqrt(-1)*(-1))+a.*b.*k0(n))+b.*exp(1) ...
  .^(sqrt(-1).*a.*c.*k0(n)).*(1+sqrt(-1).*a.*c.*k0(n))));
end
for mm = 1:length(stopPoint)
    Narray = stopPoint(mm)*NarrayPara;
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);
mask = R < Narray;
 % mask = R < Narray;
x = x(mask);
y = y(mask);
element = a*[x(:),y(:)];
origin_idx = find(element(:,1) == 0 & element(:,2) == 0, 1);

if ~isempty(origin_idx)
    element = circshift(element, -(origin_idx - 1));
end
    total_elements = size(element, 1);  % Total iterations for the inner loop
    figure(15)
    scatter(element(:,1),element(:,2))
    xlim([-a*NarrayPara, a*NarrayPara]); 
    ylim([-a*NarrayPara, a*NarrayPara]); 
    axis equal
    title('summation area')
for n = 1:length(lda)
    for ii = 2:total_elements
        Gij = tensorG(k0(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0l(:,:,mm,n) = G0l(:,:,mm,n) + Gij;
    end
end
end

ReLEFTxx = squeeze(real(G0l(1,1,:,1))) + real(G2l(:,1))*coe;
ImLEFTxx = squeeze(imag(G0l(1,1,:,1))) + imag(G2l(:,1))*coe;
ReRIGHTxx = squeeze(real(G0r(1,1,:)));
ImRIGHTxx = squeeze(imag(G0r(1,1,:)));

figure(16)
yyaxis left
plot(c,squeeze(real(G0l(1,1,:,1))))
hold on
plot(c,real(G2l(:,1))*coe,'--')
ylabel('Real part of green function')
yyaxis right
plot(c,squeeze(imag(G0l(1,1,:,1))))
hold on
plot(c,imag(G2l(:,1))*coe,'--')
ylabel('Imag part of green function')
xlabel('normalized stop summation/start integral radius')
legend('sum','int','sum','int')
title('sum vs int')
saveas(gcf, 'sum vs int.png');


figure(17)
yyaxis left
% plot(c,squeeze(real(G0l(1,1,:,1)))+real(G2l(:,1))*coe);
plot(c,ReLEFTxx/ReRIGHTxx(1))
ylabel('real part of green function')
yyaxis right
% plot(c,squeeze(imag(G0l(1,1,:,1)))+imag(G2l(:,1))*coe);
plot(c,ImLEFTxx/ImRIGHTxx(1))
ylabel('imag part of green function')
xlabel('normalized stop summation/wstart integral radius')
title('(sum + int)/Sum')
saveas(gcf, 'sum&int compare.png');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare between integral of infinite full Green function
% sum + int
NarrayPara = linspace(1e19,1e20,500);
% NarrayPara = 1e20;
b = NarrayPara;
c = 50;
G0l = zeros(3,3,length(lda));
G2lxx = zeros(length(b),length(lda));
G2lzz = zeros(length(b),length(lda));
G2l = zeros(3,3,length(b),length(lda));
for n = 1:length(lda)
G2lxx(:,n) = ((1/4).*a.^(-3).*b.^(-1).*c.^(-1).*k0(n).^(-2).*((sqrt(-1)*(-1)).*c.* ...
  exp(1).^(sqrt(-1).*a.*b.*k0(n)).*((sqrt(-1)*(-1))+a.*b.*k0(n))+b.*exp(1) ...
  .^(sqrt(-1).*a.*c.*k0(n)).*(1+sqrt(-1).*a.*c.*k0(n))));
% G2lzz(:,n) = (1/2).*a.^(-3).*c.^(-1).*k0(n).^(-2).*(sqrt(-1).*exp(1).^( ...
%   sqrt(-1).*a.*c.*k0(n)).*(sqrt(-1)+a.*c.*k0(n)));
% G2lxx(:,n) = ((1/4).*a.^(-3).*b.^(-1).*c.^(-1).*k0(n).^(-2).*(b.*exp(1) ...
%   .^(sqrt(-1).*a.*c.*k0(n)).*(1+sqrt(-1).*a.*c.*k0(n))));
G2lzz(:,n) = (1/2).*a.^(-3).*c.^(-1).*k0(n).^(-2).*(sqrt(-1).*exp(1).^( ...
  sqrt(-1).*a.*c.*k0(n)).*(sqrt(-1)+a.*c.*k0(n)));
G2l(1,1,:,n) = G2lxx(:,n)*coe;
G2l(2,2,:,n) = G2lxx(:,n)*coe;
G2l(3,3,:,n) = G2lzz(:,n)*coe;
end
Narray = c;
element = generate_hexagonal_structure(floor(Narray*2/sqrt(3)));
x = element(:,1);
y = element(:,2);
R = sqrt(x.^2 + y.^2);
mask = R < Narray;
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
        Gij = tensorG(k0(n) , 0, 0, 0, element(ii, 1), element(ii, 2), 0); % Gfarxx xx components of the tensor
        G0l(:,:,n) = G0l(:,:,n) + Gij;
    end
end

Gsum = zeros(3,3,length(b),length(lda));
for mm = 1:length(b)
Gsum(:,:,mm,:) = G0l + squeeze(G2l(:,:,mm,:));
end
% sum + int_inf

for n = 1:length(lda)
    G2infxx(n) =  ((1/4).*a.^(-3).*c.^(-1).*k0(n).^(-2).*(exp(1) ...
  .^(sqrt(-1).*a.*c.*k0(n)).*(1+sqrt(-1).*a.*c.*k0(n))));
% G2infxx(n) =  ((1/4).*a.^(-3).*c.^(-1).*k0(n).^(-2).*(exp(1) ...
%   .^(sqrt(-1).*a.*c.*k0(n)).*(1+sqrt(-1).*a.*c.*k0(n))))*coe;
G2infzz(n) = (1/2).*a.^(-3).*c.^(-1).*k0(n).^(-2).*(sqrt(-1).*exp(1).^( ...
  sqrt(-1).*a.*c.*k0(n)).*(sqrt(-1)+a.*c.*k0(n)))*coe;
end

Recomp = (real(G2infxx(length(lda))) - real(G2lxx(:,length(lda))))./real(G2infxx(length(lda)))*100;
Imcomp = (imag(G2infxx(length(lda))) - imag(G2lxx(:,length(lda))))./imag(G2infxx(length(lda)))*100;
ReGsuminfxx =  real(G2infxx') + squeeze(real(G0l(1,1,:)));
ReGsuminfzz =  real(G2infzz') + squeeze(real(G0l(3,3,:)));
ImGsuminfxx =  imag(G2infxx') + squeeze(imag(G0l(1,1,:)));
ImGsuminfzz =  imag(G2infzz') + squeeze(imag(G0l(3,3,:)));

figure(18)
yyaxis left
plot(b,Recomp)
ylabel('real part of summation green function')
yyaxis right
plot(b,Imcomp)
ylabel('imag part of summation green function')
xlabel('integral upper limit')

%%

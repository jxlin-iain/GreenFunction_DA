function cor = generate_rounded_rectangle(R, dl, ds)
xroc = [];
yroc = [];
x0 = [dl/2-R,-(dl/2-R),-(dl/2-R),dl/2-R]; % 原点横坐标
y0 = [ds/2-R,ds/2-R,-(ds/2-R),-(ds/2-R)]; % 原点纵坐标
theta_start = [0,90,180,270];   % 起始角度，单位为度
theta_end = [90,180,270,360];    % 终止角度，单位为度
n = 100;   % 线段数目
m = 50; % 直线段数目
for i = 1:length(theta_start)
theta(:,i) = linspace(theta_start(i), theta_end(i), n);   % 计算角度
end
xroc = x0 + R*cosd(theta);   % 计算x坐标
yroc = y0 + R*sind(theta);   % 计算y坐标
corroc = [reshape(xroc,1,numel(theta));reshape(yroc,1,numel(theta))];

xlength = linspace(-dl/2+R,dl/2-R,m);
ylength = linspace(-ds/2+R,ds/2-R,m);
xx = [xlength,-dl/2*ones(1,m),xlength,dl/2*ones(1,m)];
yy = [ds/2*ones(1,m),ylength,-ds/2*ones(1,m),ylength];
corlength = [xx;yy];
cor = [corroc(:,1:n),corlength(:,1:m),corroc(:,n+1:n*2),corlength(:,m+1:m*2),corroc(:,n*2+1:n*3),corlength(:,m*2+1:m*3),corroc(:,n*3+1:n*4),corlength(:,m*3+1:m*4)];

end

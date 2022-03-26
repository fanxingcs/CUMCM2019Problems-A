%% 问题一不同时刻下达到150MPa求解参数

clear,clc
global dt ftime;
ftime = 2; % 增压到150Mpa的时间
dt = 0.01; % 时间
N = 1000*ftime/dt;density = 0.850; % 密度初值
s = 0; % 计数变量
%% 寻优计算
t0s = 0:0.1:0.4; % t0的寻优列表
for t0 = t0s % 初始时高压油泵的延时
    s = s+1; % 迭代次数
    % 计算在延时为t0时的喷油时间T1和误差
    [t1(s), error(s)] = funRes(t0);  % t0：初始喷油延时，fP：目标压强
end
[Error,post] = min(error);% 选择最小的误差与对应的喷油时间
% 最优参数重新计算，记录中间的状态值
t0 = t0s(post);% 选择对应的t0与T1
T1 = t1(post);
for n = 1:N % 离散计算
    time = n * dt;  % 计算实际时间
    % 进油质量计算
    inQual = inFuel(time, density, t0, T1);
    inQ(n) = inQual; % 记录数据
    % 出油质量计算
    outQual = outFuel(time, density);
    outQ(n) = outQual;% 记录数据
    % 更新密度
    allM(n) = density*pi*5^2*500; % 计算总质量
    density = (allM(n) - outQual + inQual)/(pi*5^2*500); % 现时刻平均密度
    Denp(n) =funP2(density); % 由密度计算压强，记录值的作用
    Err1(n) = max(funP2(density)-150,0); % 计算误差
end
Err1(Err1==0)=[]; % 去除未达到的部分，消除0的影响
if isempty(Err1) 
    Err1 = inf;
end
%% 数据整理部分
figure
% subplot(3,1,3)
Err2 = sum(Err1)/N;    % 该时间内输入油量的误差
Err3 = sum(abs(Err1))/N; % 波动误差
plot((1:N)*dt,Denp,'Linewidth',1)
legend('压强变化曲线')
xlabel('时间 t(ms)')
ylabel('压强 P(Mpa)')
% set(gcf,'units','centimeters') % 标准单位：厘米
% set(gcf,'InnerPosition',[16 5 16 8]) % 在屏幕的位置，图窗的长宽
fprintf(' 延时时间为%.3f\n',t0)
fprintf('喷油角速度为%.3frad/ms\n',2*pi/T1)
fprintf(' 喷油周期为%.3f\n',T1)
fprintf(' 波动误差和为%.3f\n',Err3)
%% 计算进油工作周期
function [T1, error] = funRes(t0)
global dt ftime;
density_150 = 0.868; % 150与100Mpa对应的密度
density_100 = 0.85;
N = 1000*ftime/dt; % 总迭代次数
dep = pi*5^2*500*(density_150 - density_100)/(N*ftime); % 每次需要增加的油量
s = 0;
for Tt1 = [0.9] % 搜索法，进油工作周期
    density = 0.850; % 密度初值
    for n = 1:N % 离散计算
        time = n * dt;  % 计算实际时间
        % 进油总质量
        inQual = inFuel(time, density, t0, Tt1);
        % 出油总质量
        outQual = outFuel(time, density);
        % 更新密度
        allM = density*pi*5^2*500; % 总质量
        density = (allM - outQual + inQual - dep)/(pi*5^2*500); % 现时刻平均密度
        Err1(n) = max(funP2(density)-150,0); % 计算误差
    end
    Err1(Err1==0)=[]; % 去除未达到的部分，消除0的影响
    if isempty(Err1)
        Err1 = inf;
    end
    s = s+1;
    Tt(s) = Tt1;
    Err2(s) = sum(Err1)/N;    % 该时间内输入油量的误差
    Err3(s) = sum(abs(Err1))/N; % 波动误差
end
[~,post] = min(abs(Err2));
T1 = Tt(post);
error = Err3(post); % 波动误差
end
%% 进油计算
function inQual = inFuel(time, density, t0, Tt1)
global dt;
C = 0.85; % 流量系数
A = pi*0.7^2; % 面积
P1 = 160;
density0 = 0.871;
time = mod(time,Tt1+10); % 由于周期的存在，时间取周期的模
if time < t0
    inQual = 0;
elseif time <= t0+Tt1
    inQual = dt*C*A*sqrt(2*(P1-funP2(density))/density0)*density0;
else
    inQual = 0;
end
end
function y = funP2(density)
a0 = 0.0006492;
a2 = 1.181e-09;
a1 = -2.005e-06;
C1 = -0.217807596;
x1 = 0;
x2 = 200;
for i = 1:15
    x3 = (x1+x2)/2;
    diff = C1 + a0*x3 + a1*x3^2/2+a2*x3^3/3 - log(density);
    if diff >0
        x2 = x3;
    else
        x1 = x3;
    end
end
y = x1;
end
%% 出油量计算
function outQual = outFuel(time, density)
global dt;
time = mod(time,100); % 由于周期的存在，时间取周期的模
if time
    outQual = fPout(time)-fPout(time-dt);
    outQual = outQual * density;
else
    outQual = 0;
end
    function y = fPout(t)
        if t < 0.2
            y = 50*t^2;
        elseif t < 2.2
            y = 20*t-2;
        elseif t <=2.4
            y = 240*t-50*t^2-288+44;
        else
            y = 44;
        end
    end
end

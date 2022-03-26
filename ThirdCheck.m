%% 问题三稳定性检验

clear,clc
global dt Tc;
Tc = 1;dt = 0.1;N = 1000*Tc/dt;density = 0.850; % 密度初值
T1s = 50; % 231
% [T1, error] = funRes(T1s);  % t0：初始喷油延时，fP：目标压强
T1 = 50;
%% 重新计算最优值，记录中间结果
count = 1;
for n = 1:N % 离散计算
    time = n * dt;  % 计算实际时间
    % 每一周期过后油量增加到初始值0.5Mpa
    count0 = ceil(time/T1);
    if count0 == count
        quatyG = (8.2576-2.413)*pi*2.5^2*0.804541084;
        count0 = count;
        count = count +1;
        T1Rand = T1*(1.1-rand*0.2);
    end
    Tra(n) = T1Rand;
    % 进油 quatyG:高压油管内的燃油质量
    inQual = inFuel(time, density, quatyG, T1Rand);
    inQ(n) = inQual;
    quatyG = quatyG - inQual;
    % 出油
    outQual = outFuel(time, density);
    outQ(n) = outQual;
    % 更新密度
    allM = density*pi*5^2*500;
    density = (allM - outQual + inQual)/(pi*5^2*500);
    Denp(n) =  funP2(density);
    redu = 0;
     if Denp(n) > 100
            dP = Denp(n)-0.1;
            C = 0.85; % 流量系数
            A = pi*0.7^2;
            redu = dt*C*A*sqrt(2*dP/density)*density;
            allM = density*pi*5^2*500;
            density = (allM - redu)/(pi*5^2*500);
            Denp(n) = funP2(density);
     end
    reduQ(n) = redu;
    Err1(n) = Denp(n) - 100;
end
Err2 = sum(Err1)/N;    % 该时间内输入油量的误差
Err3 = sum(abs(Err1))/N; % 波动误差
%% 数据整理展示
% figure
% plot((1:N)*dt,Denp,'-','Linewidth',1)
% legend('压强曲线')
% xlabel('时间 t(ms)')
% ylabel('压强 P(Mpa)')
% set(gcf,'units','centimeters') % 标准单位：厘米
% set(gcf,'InnerPosition',[0 5 16 8]) % 在屏幕的位置，图窗的长宽
% figure
% plot((1:N)*dt,Err1,'-','Linewidth',1)
% legend('压强差')
% xlabel('时间 t(ms)')
% ylabel('压强差 P(Mpa)')
% set(gcf,'units','centimeters') % 标准单位：厘米
% set(gcf,'InnerPosition',[16 5 16 8]) % 在屏幕的位置，图窗的长宽
figure
set(gcf,'units','centimeters') % 标准单位：厘米
set(gcf,'InnerPosition',[0 5 13 2.9]) % 在屏幕的位置，图窗的长宽
yyaxis left
plot((1:N)*dt,inQ/dt,'b-','Linewidth',0.6)
ylim([0,60])
xlabel('时间 t(ms)')
ylabel('进油 v(mg/ms)')
yyaxis right
plot((1:N)*dt,Err1,'r-.','Linewidth',1)
ylim([-2.2,4])
xlabel('时间 t(ms)')
ylabel('压强差 P(Mpa)')
legend('进油流速','压强差','fontsize',8)

figure
set(gcf,'units','centimeters') % 标准单位：厘米
set(gcf,'InnerPosition',[0 5 13 2.9]) % 在屏幕的位置，图窗的长宽
yyaxis left
plot((1:N)*dt,outQ/dt,'b-','Linewidth',0.6)
ylim([0,25])
xlabel('时间 t(ms)')
ylabel('出油 v(mg/ms)')
yyaxis right
plot((1:N)*dt,Err1,'r-.','Linewidth',1)
ylim([-2.2,4])
xlabel('时间 t(ms)')
ylabel('压强差 P(Mpa)')
legend('出油流速','压强差','fontsize',8)

figure
set(gcf,'units','centimeters') % 标准单位：厘米
set(gcf,'InnerPosition',[0 5 13 2.9]) % 在屏幕的位置，图窗的长宽
yyaxis left
plot((1:N)*dt,reduQ/dt,'b-','Linewidth',0.6)
ylim([0,30])
xlabel('时间 t(ms)')
ylabel('阀流速 v(mg/ms)')
yyaxis right
plot((1:N)*dt,Err1,'r-.','Linewidth',1)
ylim([-2.2,4])
xlabel('时间 t(ms)')
ylabel('压强差 P(Mpa)')
legend('阀流速','压强差','fontsize',8)



fprintf('波动误差为：%.3fMpa\n',Err3)
fprintf('喷油角速度为%.3frad/ms\n',2*pi/T1)
fprintf('喷油周期为%.3fms\n',T1)
%% 计算进油工作周期
function [T1, error] = funRes(T1s)
global dt Tc;
N = 1000*Tc/dt;
s = 0;
for Tt1 = T1s % 搜索法，进油工作周期
    density = 0.850; % 密度初值
    count = 1;
    for n = 1:N % 离散计算
        time = n * dt;  % 计算实际时间
        % 每一周期过后油量增加到初始值0.5Mpa对应的油量
        count0 = ceil(time/Tt1);
        if count0 == count
            quatyG = (8.2576-2.413)*pi*2.5^2*0.804541084;
            count0 = count;
            count = count +1;
        end
        % 进油 quatyG:高压油管内的燃油质量
        inQual = inFuel(time, density, quatyG, Tt1);
        quatyG = quatyG - inQual;
        % 出油
        outQual = outFuel(time, density);
        % 更新密度
        allM = density*pi*5^2*500;
        density = (allM - outQual*2 + inQual)/(pi*5^2*500);
        thisDen = funP2(density);
        if thisDen > 100
            dP = thisDen-0.1;
            C = 0.85; % 流量系数
            A = pi*0.7^2;
            redu = dt*C*A*sqrt(2*dP/density)*density;
            allM = density*pi*5^2*500;
            density = (allM - redu)/(pi*5^2*500);
            thisDen = funP2(density);
        end
        Err1(n) = thisDen - 100;
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
function inQual = inFuel(time, density, quatyG, Tt1)
global dt;
C = 0.85; % 流量系数
A = pi*0.7^2;
time = mod(time,Tt1); % 由于周期的存在，时间取周期的模
densityThis = quatyG/(pi*2.5^2*(8.2576-funP1(time,Tt1)));
P1 = funP2(densityThis);
dP = P1-funP2(density);
if dP > 0
    inQual = dt*C*A*sqrt(2*dP/densityThis)*densityThis;
else
    inQual = 0;
end
    function y = funP1(time,Tt1)
        x = time/Tt1*2*pi;
        if x > pi
            x = 2*pi-x;
        end
        p = [0.350083971620368,-1.64948618755122,0.204612192193817,7.21756739475012];
        y =p(1)*x.^3+p(2)*x.^2+p(3)*x.^1+p(4);
    end
end
%% 出油量计算
function outQual = outFuel(time, density)
global dt;
C = 0.85;
time = mod(time,100); % 由于周期的存在，时间取周期的模
outQual = dt*C*fA(time)*sqrt(2*(funP2(density)-0.1)/density)*density;
    function y = fA(t)
        if t <= 0.33
            y = 66.71*t^3 -9.078*t^2 +0.386*t -0.00245;
        elseif t <= 2.11
            y = 1.53938040025900;
        elseif t<= 2.46
            p = [-12599.3716123347,174495.841196611,-1006433.46601700,3094198.90735766,-5347936.07171245,4926758.50356979,-1889937.68991857];
            y = p(1)*t.^6+p(2)*t.^5+p(3)*t.^4+p(4)*t.^3+p(5)*t.^2+p(6)*t.^1+p(7);
        else
            y = 0;
        end
    end
end
%% 计算一定密度下的压强
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

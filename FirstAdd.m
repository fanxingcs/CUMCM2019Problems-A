%% ����һ��ͬʱ���´ﵽ150MPa������

clear,clc
global dt ftime;
ftime = 2; % ��ѹ��150Mpa��ʱ��
dt = 0.01; % ʱ��
N = 1000*ftime/dt;density = 0.850; % �ܶȳ�ֵ
s = 0; % ��������
%% Ѱ�ż���
t0s = 0:0.1:0.4; % t0��Ѱ���б�
for t0 = t0s % ��ʼʱ��ѹ�ͱõ���ʱ
    s = s+1; % ��������
    % ��������ʱΪt0ʱ������ʱ��T1�����
    [t1(s), error(s)] = funRes(t0);  % t0����ʼ������ʱ��fP��Ŀ��ѹǿ
end
[Error,post] = min(error);% ѡ����С��������Ӧ������ʱ��
% ���Ų������¼��㣬��¼�м��״ֵ̬
t0 = t0s(post);% ѡ���Ӧ��t0��T1
T1 = t1(post);
for n = 1:N % ��ɢ����
    time = n * dt;  % ����ʵ��ʱ��
    % ������������
    inQual = inFuel(time, density, t0, T1);
    inQ(n) = inQual; % ��¼����
    % ������������
    outQual = outFuel(time, density);
    outQ(n) = outQual;% ��¼����
    % �����ܶ�
    allM(n) = density*pi*5^2*500; % ����������
    density = (allM(n) - outQual + inQual)/(pi*5^2*500); % ��ʱ��ƽ���ܶ�
    Denp(n) =funP2(density); % ���ܶȼ���ѹǿ����¼ֵ������
    Err1(n) = max(funP2(density)-150,0); % �������
end
Err1(Err1==0)=[]; % ȥ��δ�ﵽ�Ĳ��֣�����0��Ӱ��
if isempty(Err1) 
    Err1 = inf;
end
%% ����������
figure
% subplot(3,1,3)
Err2 = sum(Err1)/N;    % ��ʱ�����������������
Err3 = sum(abs(Err1))/N; % �������
plot((1:N)*dt,Denp,'Linewidth',1)
legend('ѹǿ�仯����')
xlabel('ʱ�� t(ms)')
ylabel('ѹǿ P(Mpa)')
% set(gcf,'units','centimeters') % ��׼��λ������
% set(gcf,'InnerPosition',[16 5 16 8]) % ����Ļ��λ�ã�ͼ���ĳ���
fprintf(' ��ʱʱ��Ϊ%.3f\n',t0)
fprintf('���ͽ��ٶ�Ϊ%.3frad/ms\n',2*pi/T1)
fprintf(' ��������Ϊ%.3f\n',T1)
fprintf(' ��������Ϊ%.3f\n',Err3)
%% ������͹�������
function [T1, error] = funRes(t0)
global dt ftime;
density_150 = 0.868; % 150��100Mpa��Ӧ���ܶ�
density_100 = 0.85;
N = 1000*ftime/dt; % �ܵ�������
dep = pi*5^2*500*(density_150 - density_100)/(N*ftime); % ÿ����Ҫ���ӵ�����
s = 0;
for Tt1 = [0.9] % �����������͹�������
    density = 0.850; % �ܶȳ�ֵ
    for n = 1:N % ��ɢ����
        time = n * dt;  % ����ʵ��ʱ��
        % ����������
        inQual = inFuel(time, density, t0, Tt1);
        % ����������
        outQual = outFuel(time, density);
        % �����ܶ�
        allM = density*pi*5^2*500; % ������
        density = (allM - outQual + inQual - dep)/(pi*5^2*500); % ��ʱ��ƽ���ܶ�
        Err1(n) = max(funP2(density)-150,0); % �������
    end
    Err1(Err1==0)=[]; % ȥ��δ�ﵽ�Ĳ��֣�����0��Ӱ��
    if isempty(Err1)
        Err1 = inf;
    end
    s = s+1;
    Tt(s) = Tt1;
    Err2(s) = sum(Err1)/N;    % ��ʱ�����������������
    Err3(s) = sum(abs(Err1))/N; % �������
end
[~,post] = min(abs(Err2));
T1 = Tt(post);
error = Err3(post); % �������
end
%% ���ͼ���
function inQual = inFuel(time, density, t0, Tt1)
global dt;
C = 0.85; % ����ϵ��
A = pi*0.7^2; % ���
P1 = 160;
density0 = 0.871;
time = mod(time,Tt1+10); % �������ڵĴ��ڣ�ʱ��ȡ���ڵ�ģ
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
%% ����������
function outQual = outFuel(time, density)
global dt;
time = mod(time,100); % �������ڵĴ��ڣ�ʱ��ȡ���ڵ�ģ
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

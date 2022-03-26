%% ����һ�������ڼ�����ʱ�����*_

clear,clc
global dt Tc;% ����ȫ�ֱ���dt��Tc���ֱ�Ϊʱ����ɢ���Ⱥ���󹫱�������
Tc = 2; % ��󹫱������ڣ���λs
dt = 0.01; % ʱ����ɢ����
s = 0; % ��������
%% Ѱ�ż������
t0s =  0; % Ѱ���б�
for t0 =  t0s% ��ʼʱ��ѹ�ͱõ���ʱ
    s = s+1; % ��������
    % ��������ʱΪt0ʱ������ʱ��T1�����
    [t1(s), error(s)] = funRes(t0);  % t0����ʼ������ʱ��fP��Ŀ��ѹǿ
end
% ѡ����С��������Ӧ������ʱ��
[Error,post] = min(error);% ѡ����С��������Ӧ������ʱ��
%% ���¼������Ų��� ����¼�м�״ֵ̬
t0 = t0s(post);% ��ʱʱ��
T1 = t1(post); % ��Ӧ��t1ֵ
N = 1000*Tc/dt;density = 0.850; % �ܶȳ�ֵ
for n = 1:N % ��ɢ����
    time = n * dt;  % ����ʵ��ʱ��
    % ����������
    inQual = inFuel(time, density, t0, T1);
    inQ(n) = inQual; % ��¼ֵ
    % ����������
    outQual = outFuel(time, density);
    outQ(n) = outQual;% ��¼ֵ
    % �����ܶ�
    allM(n) = density*pi*5^2*500; % ����������
    density = (allM(n) - outQual + inQual)/(pi*5^2*500); % ����ƽ���ܶ�
    Denp(n) = funP2(density); % �����ʱ��ѹ��
    Err1(n) = Denp(n) - 100; % �������
end
Err2 = sum(Err1)/N;    % ��ʱ�����������������
Err3 = sum(abs(Err1))/N; % �������
%% ��������
%����ѹ������ͼ
figure
plot((1:N)/100,Denp,'-','Linewidth',1)
legend('ѹǿ����')
xlabel('ʱ�� t(ms)')
ylabel('ѹǿ P(Mpa)')
set(gcf,'units','centimeters') % ��׼��λ������
set(gcf,'InnerPosition',[0 5 16 8]) % ����Ļ��λ�ã�ͼ���ĳ���
%�������ͼ
figure
plot((1:N)/100,Err1,'-','Linewidth',1)
legend('ѹǿ���')
xlabel('ʱ�� t(ms)')
ylabel('ѹǿ�� P(Mpa)')
set(gcf,'units','centimeters') % ��׼��λ������
set(gcf,'InnerPosition',[16 5 16 8]) % ����Ļ��λ�ã�ͼ���ĳ���
%���ƽ����������������ͼ
figure
subplot(2,1,1)
plot((1:(0.1*(N/Tc/4)))/100,inQ(1:(0.1*(N/Tc/4)))/dt,'r-','Linewidth',0.6)
xlabel('ʱ�� t(ms)')
ylabel('�����仯�� v(mg/ms)')
legend('��������������')
subplot(2,1,2)
plot((1:(0.1*(N/Tc)))/100,outQ(1:(0.1*(N/Tc)))/dt,'r-','Linewidth',0.6)
xlabel('ʱ�� t(ms)')
ylabel('�����仯�� v(mg/ms)')
legend('��������������')
fprintf('���ͽ��ٶ�Ϊ%.3frad/ms\n',2*pi/T1)
fprintf('��������Ϊ%.3fms\n',T1)
fprintf('��ʼ������ʱ%.3fms\n',t0)
fprintf('�������Ϊ��%.3fMpa\n',Err3)
%% ������͹�������
function [T1, error] = funRes(t0)
global dt Tc%����ȫ�ֱ���dt��Tc���ֱ�Ϊʱ����ɢ���Ⱥ���󹫱�������
N = 1000*Tc/dt;%һ�����ڵ���ɢʱ���ĸ���
s = 0;
for Tt1 = 27:30 % �����������͹�������
    density = 0.850; % �ܶȳ�ֵ
    for n = 1:N % ��ɢ����
        time = n * dt;  % ����ʵ��ʱ��
        % ����
        inQual = inFuel(time, density, t0, Tt1*dt);%InFuel�ǽ��ͼ��㺯��
        % ����
        outQual = outFuel(time, density);%outFuel�ǳ��ͼ��㺯��
        % �Ƿ�����ѹ
        
        % �����ܶ�
        allM = density*pi*5^2*500;% ����������
        density = (allM - outQual + inQual)/(pi*5^2*500);% ����ƽ���ܶ�
        Err1(n) = funP2(density) - 100;% �������
    end
    s = s+1;
    Tt(s) = Tt1*dt;
    Err2(s) = sum(Err1)/N;    % ��ʱ�����������������
    Err3(s) = sum(abs(Err1))/N; % �������
end
[~,post] = min(abs(Err2));% ��С���
T1 = Tt(post);
error = Err3(post); % �������
end
%% ���ͼ���
function inQual = inFuel(time, density, t0, Tt1)
global dt;% ����ȫ�ֱ���dt��Ϊʱ����ɢ����
C = 0.85; % ����ϵ��
A = pi*0.7^2;% A�ڵ����
P1 = 160;% ��ѹ�ͱ���A���ṩ�ĺ�ѹ
density0 = 0.871;% ��ʼ�ܶ�
time = mod(time,Tt1+10); % �������ڵĴ��ڣ�ʱ��ȡ���ڵ�ģ

if time < t0
    inQual = 0;
elseif time <= t0+Tt1% һ����������
    inQual = dt*C*A*sqrt(2*(P1-funP2(density))/density0)*density0;
else
    inQual = 0;
end
end
%% �ܶ�ѹǿת������
function y = funP2(density)
a0 = 0.0006492;
a2 = 1.181e-09;
a1 = -2.005e-06;
C1 = -0.217807596;
x1 = 0;
x2 = 200;
for i = 1:15
    x3 = (x1+x2)/2;
    diff = C1 + a0*x3 + a1*x3^2/2+a2*x3^3/3 - log(density);%�ɸ���������϶���
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
global dt;% ����ȫ�ֱ���dt��Ϊʱ����ɢ����
time = mod(time,100); % �������ڵĴ��ڣ�ʱ��ȡ���ڵ�ģ
if time
    outQual = fPout(time)-fPout(time-dt);
    outQual = outQual * density;
else
    outQual = 0;
end
    function y = fPout(t)
        if t<50-0.4
            y = 0;
        elseif t < 0.2+50-0.4
            t = t-50+0.4;
            y = 50*t^2;
        elseif t < 2.2+50-0.4
            t = t-50+0.4;
            y = 20*t-2;
        elseif t <=2.4+50-0.4
            t = t-50+0.4;
            y = 240*t-50*t^2-288+44;
        else
            y = 44;
        end
    end
end

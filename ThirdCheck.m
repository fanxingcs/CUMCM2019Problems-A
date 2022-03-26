%% �������ȶ��Լ���

clear,clc
global dt Tc;
Tc = 1;dt = 0.1;N = 1000*Tc/dt;density = 0.850; % �ܶȳ�ֵ
T1s = 50; % 231
% [T1, error] = funRes(T1s);  % t0����ʼ������ʱ��fP��Ŀ��ѹǿ
T1 = 50;
%% ���¼�������ֵ����¼�м���
count = 1;
for n = 1:N % ��ɢ����
    time = n * dt;  % ����ʵ��ʱ��
    % ÿһ���ڹ����������ӵ���ʼֵ0.5Mpa
    count0 = ceil(time/T1);
    if count0 == count
        quatyG = (8.2576-2.413)*pi*2.5^2*0.804541084;
        count0 = count;
        count = count +1;
        T1Rand = T1*(1.1-rand*0.2);
    end
    Tra(n) = T1Rand;
    % ���� quatyG:��ѹ�͹��ڵ�ȼ������
    inQual = inFuel(time, density, quatyG, T1Rand);
    inQ(n) = inQual;
    quatyG = quatyG - inQual;
    % ����
    outQual = outFuel(time, density);
    outQ(n) = outQual;
    % �����ܶ�
    allM = density*pi*5^2*500;
    density = (allM - outQual + inQual)/(pi*5^2*500);
    Denp(n) =  funP2(density);
    redu = 0;
     if Denp(n) > 100
            dP = Denp(n)-0.1;
            C = 0.85; % ����ϵ��
            A = pi*0.7^2;
            redu = dt*C*A*sqrt(2*dP/density)*density;
            allM = density*pi*5^2*500;
            density = (allM - redu)/(pi*5^2*500);
            Denp(n) = funP2(density);
     end
    reduQ(n) = redu;
    Err1(n) = Denp(n) - 100;
end
Err2 = sum(Err1)/N;    % ��ʱ�����������������
Err3 = sum(abs(Err1))/N; % �������
%% ��������չʾ
% figure
% plot((1:N)*dt,Denp,'-','Linewidth',1)
% legend('ѹǿ����')
% xlabel('ʱ�� t(ms)')
% ylabel('ѹǿ P(Mpa)')
% set(gcf,'units','centimeters') % ��׼��λ������
% set(gcf,'InnerPosition',[0 5 16 8]) % ����Ļ��λ�ã�ͼ���ĳ���
% figure
% plot((1:N)*dt,Err1,'-','Linewidth',1)
% legend('ѹǿ��')
% xlabel('ʱ�� t(ms)')
% ylabel('ѹǿ�� P(Mpa)')
% set(gcf,'units','centimeters') % ��׼��λ������
% set(gcf,'InnerPosition',[16 5 16 8]) % ����Ļ��λ�ã�ͼ���ĳ���
figure
set(gcf,'units','centimeters') % ��׼��λ������
set(gcf,'InnerPosition',[0 5 13 2.9]) % ����Ļ��λ�ã�ͼ���ĳ���
yyaxis left
plot((1:N)*dt,inQ/dt,'b-','Linewidth',0.6)
ylim([0,60])
xlabel('ʱ�� t(ms)')
ylabel('���� v(mg/ms)')
yyaxis right
plot((1:N)*dt,Err1,'r-.','Linewidth',1)
ylim([-2.2,4])
xlabel('ʱ�� t(ms)')
ylabel('ѹǿ�� P(Mpa)')
legend('��������','ѹǿ��','fontsize',8)

figure
set(gcf,'units','centimeters') % ��׼��λ������
set(gcf,'InnerPosition',[0 5 13 2.9]) % ����Ļ��λ�ã�ͼ���ĳ���
yyaxis left
plot((1:N)*dt,outQ/dt,'b-','Linewidth',0.6)
ylim([0,25])
xlabel('ʱ�� t(ms)')
ylabel('���� v(mg/ms)')
yyaxis right
plot((1:N)*dt,Err1,'r-.','Linewidth',1)
ylim([-2.2,4])
xlabel('ʱ�� t(ms)')
ylabel('ѹǿ�� P(Mpa)')
legend('��������','ѹǿ��','fontsize',8)

figure
set(gcf,'units','centimeters') % ��׼��λ������
set(gcf,'InnerPosition',[0 5 13 2.9]) % ����Ļ��λ�ã�ͼ���ĳ���
yyaxis left
plot((1:N)*dt,reduQ/dt,'b-','Linewidth',0.6)
ylim([0,30])
xlabel('ʱ�� t(ms)')
ylabel('������ v(mg/ms)')
yyaxis right
plot((1:N)*dt,Err1,'r-.','Linewidth',1)
ylim([-2.2,4])
xlabel('ʱ�� t(ms)')
ylabel('ѹǿ�� P(Mpa)')
legend('������','ѹǿ��','fontsize',8)



fprintf('�������Ϊ��%.3fMpa\n',Err3)
fprintf('���ͽ��ٶ�Ϊ%.3frad/ms\n',2*pi/T1)
fprintf('��������Ϊ%.3fms\n',T1)
%% ������͹�������
function [T1, error] = funRes(T1s)
global dt Tc;
N = 1000*Tc/dt;
s = 0;
for Tt1 = T1s % �����������͹�������
    density = 0.850; % �ܶȳ�ֵ
    count = 1;
    for n = 1:N % ��ɢ����
        time = n * dt;  % ����ʵ��ʱ��
        % ÿһ���ڹ����������ӵ���ʼֵ0.5Mpa��Ӧ������
        count0 = ceil(time/Tt1);
        if count0 == count
            quatyG = (8.2576-2.413)*pi*2.5^2*0.804541084;
            count0 = count;
            count = count +1;
        end
        % ���� quatyG:��ѹ�͹��ڵ�ȼ������
        inQual = inFuel(time, density, quatyG, Tt1);
        quatyG = quatyG - inQual;
        % ����
        outQual = outFuel(time, density);
        % �����ܶ�
        allM = density*pi*5^2*500;
        density = (allM - outQual*2 + inQual)/(pi*5^2*500);
        thisDen = funP2(density);
        if thisDen > 100
            dP = thisDen-0.1;
            C = 0.85; % ����ϵ��
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
    Err2(s) = sum(Err1)/N;    % ��ʱ�����������������
    Err3(s) = sum(abs(Err1))/N; % �������
end
[~,post] = min(abs(Err2));
T1 = Tt(post);
error = Err3(post); % �������
end
%% ���ͼ���
function inQual = inFuel(time, density, quatyG, Tt1)
global dt;
C = 0.85; % ����ϵ��
A = pi*0.7^2;
time = mod(time,Tt1); % �������ڵĴ��ڣ�ʱ��ȡ���ڵ�ģ
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
%% ����������
function outQual = outFuel(time, density)
global dt;
C = 0.85;
time = mod(time,100); % �������ڵĴ��ڣ�ʱ��ȡ���ڵ�ģ
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
%% ����һ���ܶ��µ�ѹǿ
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

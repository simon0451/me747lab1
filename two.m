clear all;
close all;

%% 2.1 b) Finding Experimental damping ratio and undamped natural frequency
filename1 = '2.1.4 Data.lvm';
output1 = importdata(filename1,'\t',33);
responsedata1 = output1.data(:,3:4);

% find peaks
[maxtab,mintab] = peakdet(responsedata1(:,2),0.05);

% check peaks
figure(1);
grid on;
hold on;
plot(responsedata1(:,1),responsedata1(:,2));
plot(responsedata1(maxtab(:,1),1),maxtab(:,2),'O');
title('Experimental Step Response Data','FontSize',14);
xlabel('Time (s)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','Data','Peaks');

% Log decrement method
n = 2; % next peak is 5 periods away
sigma = (1/(n-1))*log(maxtab(1,2)/maxtab(2,2));
zeta = sigma/sqrt(4*pi^2+(sigma)^2); % the damping ratio
T = diff(responsedata1(maxtab(:,1),1));
T = T(1); % The first damped period
wd = 2*pi/T; % The damped frequency in rad/s -> 2pi rad/cyc div s/cyc
wn = wd/(sqrt(1-zeta^2)); % Natural frequency in rad/s
wnHz = wn/(2*pi); % Natural frequency in Hz -> rad/s div 2pi rad/cyc

%% c)
% From derived differential equation, LC = 1/wn^2, L/R = 2*zeta/wn
R = 20200; % R = 20.2k ohms measured in lab
L = R*2*zeta/wn; % result = 1.3 mH
C = 1/wn^2/L; % result = 179.21 muF

sys = tf(1,[L*C,L/R,1]);

%% d) Step Responses
% creates the step offset conditions
opt1 = stepDataOptions('InputOffset',-2,'StepAmplitude',4);

% step response
[V1,t1] = step(sys,opt1);

figure(2);
hold on;
grid on;
xlim([-0.0025 0.035]);
ylim([-2 6]);
plot(t1,V1,'--');
plot(responsedata1(:,1),responsedata1(:,2));
title('Simulated Step Response of the Second Order Circuit','FontSize',14);
xlabel('Time (ms)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
legend('Location','best','Simulated Response','Experimental Response');

%% 2.2 a) - done in excel sheet

%% b)
omega_i = [10,16,80,158,317,476,500,634,1000,1268,5000,7000]; % Hz
dB = [0,0.043647778,0.4477155,2.209589084,8.93064944,-2.579369816,... % dB
    -3.777147669,-9.11551413,-17.79041092,-22.12060235,-46.40073812,-53.18542508];
phi = [-0.2,-0.3,-1.5,-5.9,-89,-162,-164,-169,-174,-175,-178,-178]; % deg

% retrieve bode plot data from bode(sys)
[mag,phase,wout] = bode(sys);
phase = squeeze(phase); % reduce the 3D matrix into 2D
mag = squeeze(mag); % reduce the 3D matrix into 2D
wout = wout/(2*pi); % convert from rad/s to Hz

% convert mag1 to decibels
magdb = 20*log10(mag);

figure(3);
subplot(2,1,1);
semilogx(wout,magdb,'kO');
hold on;
grid on;
semilogx(omega_i,dB,'b');
xlim([10 10000]);
ylim([-40 25]);
ylabel('Amplitude Ratio (dB)','FontSize',12);
title('Simulated and Experimental Bode Response Plot','FontSize',14);
legend('Location','best','Simulated Bode Plot','Experimental Bode Plot');

subplot(2,1,2);
semilogx(wout,phase,'kO');
hold on;
grid on;
semilogx(omega_i,phi,'b');
xlim([10 10000]);
ylim([-180 0]);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Phase Shift (deg)','FontSize',12);
legend('Location','best','Simulated Bode Plot','Experimental Bode Plot');

%% c)
[maxPeak,maxInd] = max(magdb); % maximum point on the measured bode plot
wnFromBode = wout(maxInd)*2*pi; % frequency where maxPeak happened
zetaExp = 1/(2*maxPeak); % Approximate height = 1/2zeta

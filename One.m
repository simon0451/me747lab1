clear all;
close all;

%% Reading lvm data for first output
filename1 = '1.1.5 Data.lvm';
output1 = importdata(filename1,'\t',33);
responsedata1 = output1.data(:,3:4);

% Output 1 from 1.1.5
% Observed from data, t = -3.13e-7 is the threshold
% Takes first 5 values to find the slope.
% Second method

threshold1 = -3.13e-7; % this is the first observable starting point of the response curve
responseFinalValue1 = mean(responsedata1(4950:end,2));
responseInitialValue1 = mean(responsedata1(1:350,2));

[onetauX1,onetauX2,onetauY1,onetauY2,onetauPer1,onetauPer2,p1] = ...
    FindTaus(responsedata1,threshold1,responseFinalValue1,responseInitialValue1);

onetauTime1 = onetauX1 - threshold1;
onetauTime2 = onetauX2 - threshold1;

% The initial slope line
t1 = linspace(0,onetauTime1,100);
o1 = t1*p1(1) + p1(2);

% Plot the results
figure;
hold on;
grid on;
xlim([-1 5]);
ylim([-2.5 2.5]);
% plot the original data with initial time offset adjusted
plot((responsedata1(:,1)-threshold1)*1e3,responsedata1(:,2),'-.');
% plot the line
plot(t1*1e3,o1);
% plot the tau location on the curve from the initial slope method
plot(onetauTime1*1e3,onetauY1,'O');
% plot the tau location on the curve from the second method
plot(onetauTime2*1e3,onetauY2,'X','markers',12);
% plot a drop down from the tip of the line to the tau location
plot([onetauTime1*1e3 onetauTime1*1e3],[onetauY1 o1(end)],'--');
legend('Location','best','Original response data','Initial slope line',...
    'Tau point from initial slope method','Tau point from 63.2% method',...
    'Drop down indicator line for initial slope method');
title('First Order Step Response 0 Offset Analysis','FontSize',14);
xlabel('Time (ms)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);

%% Reading lvm data for second output
filename2 = '1.1.6 Data.lvm';
output2 = importdata(filename2,'\t',33);
responsedata2 = output2.data(:,3:4);

% Output 1 from 1.1.6
% observed from data, t = 2.717e-7 is the threshold
% Takes first 5 values to find the slope.

threshold2 = 2.717e-7; % this is the first value higher than the initial noise
responseFinalValue2 = mean(responsedata2(4950:end,2));
responseInitialValue2 = mean(responsedata2(1:350,2));

[twotauX1,twotauX2,twotauY1,twotauY2,twotauPer1,twotauPer2,p2] = ...
    FindTaus(responsedata2,threshold2,responseFinalValue2,responseInitialValue2);

twotauTime1 = twotauX1 - threshold2;
twotauTime2 = twotauX2 - threshold2;

% The initial slope line
t2 = linspace(0,twotauX1,100);
o2 = t2*p2(1) + p2(2);

% Plot the results
figure;
hold on;
grid on;
xlim([-1 5]);
ylim([-3.5 5.5]);
% plot the original data with initial time offset adjusted
plot((responsedata2(:,1)-threshold2)*1e3,responsedata2(:,2),'-.');
% plot the line
plot(t2*1e3,o2);
% plot the tau location on the curve from the initial slope method
plot(twotauTime1*1e3,twotauY1,'O');
% plot the tau location on the curve from the second method
plot(twotauTime2*1e3,twotauY2,'X','markers',12);
% plot a drop down from the tip of the line to the tau location
plot([twotauTime1*1e3 twotauTime1*1e3],[twotauY1 o2(end)],'--');
legend('Location','best','Original response data','Initial slope line',...
    'Tau point from initial slope method','Tau point from 63.2% method',...
    'Drop down indicator line for initial slope method');
title('First Order Step Response 1V Offset Analysis','FontSize',14);
xlabel('Time (ms)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);

%% Theoretical Simulation
R = 11.1e3;
C = 58.77e-9;
sys = tf(1,[R*C 1]);

opt = stepDataOptions('InputOffset',-2,'StepAmplitude',4);

[V,t] = step(sys,opt);

% Find tau
for i = 1:length(V)
    if (V(i) >= 0.632*(V(end)-V(1))+V(1))
        tauX = t(i);
        tauY = V(i);
        break;
    end
end

figure;
hold on;
grid on;
xlim([-1 5]);
ylim([-2.5 2.5]);
plot(t*1e3,V);
plot(tauX*1e3,tauY,'*','markers',14);
plot([-1 5],[tauY tauY],'k--');
title('Simulated Response of the First Order Circuit','FontSize',14);
xlabel('Time (ms)','FontSize',12);
ylabel('Voltage (V)','FontSize',12);
text(2,-0.5,['Tau = ',num2str(tauX*1e3),' ms']);
text(2,-0.8,['Tau Voltage = ',num2str(tauY),' V']);

figure;
hold on;
grid on;
xlim([-1 5]);
ylim([-2.5 2.5]);
plot(t*1e3,V);
plot((responsedata1(:,1)-threshold1)*1e3,responsedata1(:,2),'-.');
legend('Location','best','Theoretical Simulation','Experimental Result');
title('Comparing Theoretical and Experimental Step Response','FontSize',14);
xlabel('Time (ms)','FontSize',12);
ylabel('VOltage (V)','FontSize',12);
clear all;
close all;

%% Reading lvm data for first output
filename1 = '1.1.5 Data.lvm';
output1 = importdata(filename1,'\t',33);
responsedata1 = output1.data(:,3:4);

% Output 1 from 1.1.5
% observed from data, y = -1.949 is the threshold
% observed from data, between threshold and y = -1.471 is a good segment
% to take the initial slope.
% Chocies are: observe a segment or use a set number of data points

threshold1 = -3.13e-7; % this is the first observable starting point of the response curve
responseFinalValue1 = mean(responsedata1(4950:end,2));
responseInitialValue1 = mean(responsedata1(1:350,2));

[onetauX1,onetauX2,onetauY1,onetauY2,onetauPer1,onetauPer2] = ...
    FindTaus(responsedata1,threshold1,responseFinalValue1,responseInitialValue1);

onetauTime1 = onetauX1 - threshold1;
onetauTime2 = onetauX2 - threshold1;

%% Reading lvm data for second output
filename2 = '1.1.6 Data.lvm';
output2 = importdata(filename2,'\t',33);
responsedata2 = output2.data(:,3:4);

% Output 1 from 1.1.6
% observed from data, y = -2.91 is the threshold
% observed from data, between threshold and y = -2.721 is a good segment
% to take the initial slope.
% Chocies are: observe a segment or use a set number of data points

threshold2 = 2.717e-7; % this is the first value higher than the initial noise
responseFinalValue2 = mean(responsedata2(4950:end,2));
responseInitialValue2 = mean(responsedata2(1:350,2));

[twotauX1,twotauX2,twotauY1,twotauY2,twotauPer1,twotauPer2] = ...
    FindTaus(responsedata2,threshold2,responseFinalValue2,responseInitialValue2);

twotauTime1 = twotauX1 - threshold2;
twotauTime2 = twotauX2 - threshold2;
clear all;
close all;

%% Reading lvm data for first output
filename1 = '2.1.4 Data.lvm';
output1 = importdata(filename1,'\t',33);
responsedata1 = output1.data(:,3:4);

plot(responsedata1(:,1),responsedata1(:,2));


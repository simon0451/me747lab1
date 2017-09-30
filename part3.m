%This script is responsible for Part 3 of the analysis
%Author: Simon Popecki
%Date: 30 September 2017
%Class: ME 747

%Part 3 concerns the response of the first and second order systems to a
%random noise input. The experimental values are compared to the
%theoretical predictions.

%% Data Importation
clear all;
close all;

load Part3Data.mat
%DATA FORMAT:
%For amplitude response plots, the X-value from the raw data is frequency (Hz)
%The Y-value is the corresponding magnitude (dB)
%For phase response plots, the X-vlue from the raw data is frequency (Hz)
%The Y-value is the corresponding phase angle (degrees)

%Signal Conditioning
%The raw bode plots from LabView have a lot of noise past 10 kHz (1st
%order) and 1.4 kHz (2nd order)

%Trimming the 1st order bode plots at ~10 kHz
frequencyAR1 = frequencyAR1(1:1000);
magnitudeAR1 = magnitudeAR1(1:1000);
frequencyPS1 = frequencyPS1(1:1000);
phaseAnglePS1 = phaseAnglePS1(1:1000);



%% Freqeuncy Response of a First Order System Using LabView 

figure(1)
subplot(2,1,1)
semilogx(frequencyAR1,magnitudeAR1)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
subplot(2,1,2)
semilogx(frequencyPS1,phaseAnglePS1)
xlabel('Frequency (Hz)')
ylabel('Phase Angle (°)')
grid on


%% Frequency Response of a Second Order System Using LabView
figure(2)
subplot(2,1,1)
semilogx(frequencyAR2,magnitudeAR2)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on
subplot(2,1,2)
semilogx(frequencyPS2,phaseAnglePS2)
xlabel('Frequency (Hz)')
ylabel('Phase Angle (°)')
grid on















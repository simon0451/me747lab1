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

%Trimming the 1st order bode plot at ~10 kHz
cutoff1 = 1000; %frequency index cutoff number, first order
frequencyAR1 = frequencyAR1(1:cutoff1);
magnitudeAR1 = magnitudeAR1(1:cutoff1);
frequencyPS1 = frequencyPS1(1:cutoff1);
phaseAnglePS1 = phaseAnglePS1(1:cutoff1);

%Trimming the 2nd order bode plot at 1.4 kHz
cutoff2 = 500; %frequency index cutoff number, second order
frequencyAR2 = frequencyAR2(1:cutoff2);
magnitudeAR2 = magnitudeAR2(1:cutoff2);
frequencyPS2 = frequencyPS2(1:cutoff2);
phaseAnglePS2 = phaseAnglePS2(1:cutoff2);


%% Freqeuncy Response of a First Order System Using LabView 
%First order known values
R = 11.1e3; %Ohms, the resistance of the resistor
C = 58.77e-9; %Farrads, the capacitance of the capacitor
sys1 = tf(1,[R*C 1]); %The transfer function of the 1st order system

[mag,theoreticalPhaseAngle1,theoreticalFrequency1] = bode(sys1); %Extracting data points from system object, sys1
theoreticalPhaseAngle1 = squeeze(theoreticalPhaseAngle1); %force array
mag1 = squeeze(mag); %force array
TheoreticalMagnitudeAR1 = 20*log10(mag1); %converting from absolute magnitude to decibels
theoreticalFrequency1 = theoreticalFrequency1/(2*pi); %convert from rad/s to Hz



figure(1)
subplot(2,1,1)
semilogx(frequencyAR1,magnitudeAR1,theoreticalFrequency1,TheoreticalMagnitudeAR1,'o')
title('Theoretical and Experimental Frequency Response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Experimental Values','Theoretical Values','location','southwest')
grid on
subplot(2,1,2)
semilogx(frequencyPS1,phaseAnglePS1,theoreticalFrequency1,theoreticalPhaseAngle1,'o')
xlabel('Frequency (Hz)')
ylabel('Phase Angle (°)')
legend('Experimental Values','Theoretical Values','location','southwest')
grid on


%% Frequency Response of a Second Order System Using LabView
figure(2)
subplot(2,1,1)
semilogx(frequencyAR2,magnitudeAR2)
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Experimental Values','location','southwest')
grid on
subplot(2,1,2)
semilogx(frequencyPS2,phaseAnglePS2)
xlabel('Frequency (Hz)')
ylabel('Phase Angle (°)')
legend('Experimental Values','location','southwest')
grid on















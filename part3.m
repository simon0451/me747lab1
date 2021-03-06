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
%First order known values (R and C get overwritten)
R = 11.1e3; %Ohms, the resistance of the resistor
C = 58.77e-9; %Farrads, the capacitance of the capacitor
sys1 = tf(1,[R*C 1]); %The transfer function of the 1st order system

[mag,theoreticalPhaseAngle1,theoreticalFrequency1] = bode(sys1); %Extracting data points from system object, sys1
theoreticalPhaseAngle1 = squeeze(theoreticalPhaseAngle1); %force array
mag1 = squeeze(mag); %force array
TheoreticalMagnitudeAR1 = 20*log10(mag1); %converting from absolute magnitude to decibels
theoreticalFrequency1 = theoreticalFrequency1/(2*pi); %convert from rad/s to Hz

%1/Tau (the time constant) occurs at the point where phase is 45 degrees
%omega is ~250 Hz
locationMap = find(phaseAnglePS1<=-45); %finding the indices of phase angles that meet the condition <= -45 degrees
freqInd = locationMap(1); %finding the index of the first value meeting the <= -45 degrees condition
omega1 = frequencyPS1(freqInd); %Hz, 1/tau - the natural frequency of the first order system
tau1 = 1/(omega1*6.28); %seconds, the time constant of the first order system (multiplied by 6.28 to convert to radians per second
disp('Circuit 1 time constant (seconds):')
disp(tau1)
C1 = 1/(R*omega1*6.28); %Farrads, the calculated capacitance of the first order system using data from LabView
disp('Circuit 1 capacitance (Farrads):')
disp(C1)

figure(1)
subplot(2,1,1)
semilogx(frequencyAR1,magnitudeAR1,theoreticalFrequency1,TheoreticalMagnitudeAR1,'o')
% title('Theoretical and Experimental Frequency Response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Experimental Values','Theoretical Values','location','southwest')
grid on
xmin = 0;
xmax = 100000;
ymin = -40;
ymax = 5;
axis ([xmin xmax ymin ymax])
subplot(2,1,2)
semilogx(frequencyPS1,phaseAnglePS1,theoreticalFrequency1,theoreticalPhaseAngle1,'o')
xlabel('Frequency (Hz)')
ylabel('Phase Angle (�)')
legend('Experimental Values','Theoretical Values','location','southwest')
grid on
xmin = 0;
xmax = 100000;
ymin = -100;
ymax = 10;
axis ([xmin xmax ymin ymax])


%% Frequency Response of a Second Order System Using LabView
%Second order known values
L = 1.272690763057150; %H, inductance
R = 20200; %Ohms, the resistance of the resistor
C = 2.030234051171926e-07; %Farrads, the capacitance of the capacitor
sys2 = tf(1,[L*C L/R 1]); %The transfer function of the 1st order system

[mag,theoreticalPhaseAngle2,theoreticalFrequency2] = bode(sys2); %Extracting data points from system object, sys1
theoreticalPhaseAngle2 = squeeze(theoreticalPhaseAngle2); %force array
mag2 = squeeze(mag); %force array
TheoreticalMagnitudeAR2 = 20*log10(mag2); %converting from absolute magnitude to decibels
theoreticalFrequency2 = theoreticalFrequency2/(2*pi); %convert from rad/s to Hz

%Finding Experimental Damping Ratio
peakdB = max(magnitudeAR2); %finding the magnitude of the peak
zeta = 1/(peakdB*2); %calculating the damping ratio given the peak magnitude
disp('Circuit 2 damping ratio:')
disp(zeta)

%Finding Experimental Natural Frequency
[~,indexMaxMag] = max(magnitudeAR2); %finding the index of the max amplitude
natFrequency = frequencyAR2(indexMaxMag); %Hz, finding the frequency at the break point
disp('Circuit 2 natural frequency (Hz):')
disp(natFrequency)

%Calculating L and C values from LabView experimental results
LC = 1/((natFrequency*6.28)^2);
L = R*(2*zeta/(natFrequency*6.28)); %H, calculating inductance
C = LC/L;
disp('Circuit 2 Inductance (H):')
disp(L)
disp('Circuit 2 Resistance (Ohms):')
disp(R)
disp('Circuit 2 Capacitance (F):')
disp(C)

figure(2)
subplot(2,1,1)
semilogx(frequencyAR2,magnitudeAR2,theoreticalFrequency2,TheoreticalMagnitudeAR2,'o')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
legend('Experimental Values','Theoretical Values','location','southwest')
grid on
xmin = 0;
xmax = 10000;
ymin = -60;
ymax = 25;
axis ([xmin xmax ymin ymax])
subplot(2,1,2)
semilogx(frequencyPS2,phaseAnglePS2,theoreticalFrequency2,theoreticalPhaseAngle2,'o')
xlabel('Frequency (Hz)')
ylabel('Phase Angle (�)')
legend('Experimental Values','Theoretical Values','location','southwest')
grid on
xmin = 0;
xmax = 10000;
ymin = -200;
ymax = 10;
axis ([xmin xmax ymin ymax])














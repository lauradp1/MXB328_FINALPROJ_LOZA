clear; clc; close all
%% Import Data
rainfall = readtable('MackayAlert.csv');
%% Climate Modelling Plots
[cosineRain, fourierRain, r_f, average] = RainfallModelling(rainfall);
% Average Rainfall Plot
figure;
plot(average*10^(-3),'LineWidth',2)
hold on 
fplot(cosineRain, [0 366],'LineWidth',2)
title('Rainfall Averaged from 2012 to 2020 Compared with Cosine Approximation','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
ylabel('Rainfall (m)','FontSize',20,'Interpreter','LaTeX')
xlim([0 366])
ylim([0 max(average*10^(-3))])
lgd = legend('Averaged Plot','Cosine Approximation');
lgd.Interpreter = 'latex';
lgd.FontSize = 16;
%% Average Plot with Fourier Approximation
t = linspace(1,366,366);
figure;
plot(t,average,'LineWidth',2)
hold on
fplot(fourierRain, [0 366],'r--','LineWidth',2)
title('Rainfall Averaged from 2012 to 2020','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
xlim([0 366])
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')
lgd = legend('Averaged Data','Fourier Approximation');
lgd.Interpreter = 'LaTeX';
lgd.FontSize = 16;
%% Advanced Fourier
clear; clc; close all
%% Import Data
rainfall = readtable('MackayAlert.csv');
%% Climate Modelling Plots
[cosineRain, fourierRain, r_f, average, fourierFlood, cosineFlood, RainfallModel] = RainfallModelling(rainfall);
%% Average Rainfall Plot
figure;
plot(average,'LineWidth',2)
hold on 
fplot(cosineRain, [0 366],'LineWidth',2)
hold on
yline(r_f,'k-','LineWidth',2)
title('Rainfall Averaged from 2012 to 2020 Compared with Cosine and Constant Approximations','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')
xlim([0 366])
ylim([0 (max(average)+1)])
lgd = legend('Averaged Plot','Cosine Approximation','Constant Rainfall');
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


%% Cosine Flood
figure
plot(t, RainfallModel(:,1),'b','LineWidth',2)
hold on
fplot(cosineFlood, [0 366],'r--','LineWidth',2) 
title('Cosine Approximation Compared with Flood Year Data','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
xlim([0 366])
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')
lgd = legend('2011 Data','Cosine Approximation');
lgd.Interpreter = 'latex';

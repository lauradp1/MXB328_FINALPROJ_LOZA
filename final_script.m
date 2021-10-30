clear; clc; close all
%% Import Data
rainfall = readtable('MackayAlert.csv');
%% Climate Modelling Plots
[years, cosineRain, fourierRain, r_f, average, fourierFlood, cosineFlood, RainfallModel, steadystateFlood, steadystateAveraged, FourierAverage, FourierFlood] = RainfallModelling(rainfall);
%% Initial Yearly Rainfall Plot
figure;
t = linspace(1,366,366);
plot(t,RainfallModel(:,1),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,2),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,3),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,4),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,5),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,6),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,7),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,8),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,9),'LineWidth',1.5)
hold on
plot(t,RainfallModel(:,10),'LineWidth',1.5)
xlim([0 366])
title('Rainfall from 2011 to 2020 Compared','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')
lgd = legend('2011','2012','2013','2014','2015','2016','2017','2018','2019','2020');
lgd.FontSize = 16;
lgd.Interpreter = 'latex';
lgd.Location = 'eastoutside';
%% How do I make this a loop ???
% t = linspace(1,366,366);
% figure
% for i = 1:length(years)
%     plot(t,RainfallModel(:,i))
%     hold on
% end
% xlim([0 366])
%% Average Rainfall Plot
figure;
plot(average,'b','LineWidth',2)
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
plot(t,average,'b','LineWidth',2)
hold on
plot(t,FourierAverage,'r--','LineWidth',2)
xlim([0 366])
title('Rainfall Averaged from 2012 to 2020','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')
lgd = legend('Averaged Data','Fourier Approximation');
lgd.Interpreter = 'LaTeX';
lgd.FontSize = 16;
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
lgd.FontSize = 16;
%% Fourier Flood 
figure
plot(t, RainfallModel(:,1),'b','LineWidth',2)
hold on
plot(t,FourierFlood,'r--','LineWidth',2) 
title('Fourier Approximation Compared with Flood Year Data','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
xlim([0 366])
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')
lgd = legend('2011 Data','Fourier Approximation');
lgd.Interpreter = 'latex';
lgd.FontSize = 16;

close all;
clc;
%clear all;
%% Import Data
rainfall = readtable('MackayAlert.csv');
%% Remove unneccesary columns and convert to matrix
rainfall(:,[1,2,7,8]) = [];
rainfall = table2array(rainfall);
%% Remove all other years outside of 2012-2020
rainfall = rainfall(rainfall(:,1)<2021 & rainfall(:,1)>=2012,:);
%% Change format of data
years = unique(rainfall(:,1));
rainfallData = zeros(366,length(years));
for Average = 1:length(years)
    yearData = rainfall(rainfall(:,1)==years(Average),4);
    if length(yearData) ~= 366
        yearData = [yearData(1:59);NaN;yearData(60:end)];
    end
    rainfallData(:,Average) = yearData;
end
%% Statistics
meanRainfall = zeros(length(years),1);
for i = 1:length(years)
    meanRainfall(i) = mean(rainfallData(~isnan(rainfallData(:,i)),i));
end
maxRainfall = max(rainfallData);
minRainfall = min(rainfallData);
%% Overall
MAX = max(maxRainfall);
MIN = min(minRainfall);
MEAN = mean(meanRainfall);
%% Cosine Approximation
x1 = (1:366);
a1 = 0.5*MAX;
b1 = 2*pi/365;
d1 = a1;
y1 = a1*cos(b1*x1)+d1;
%% Piecewise Approximation
% Cosine
x2 = (1:182);
a2 = 0.5*MAX;
b2 = 2*pi/365;
d2 = a2;
y2 = a2*cos(b2*x2)+d2;
% Sine 
x3 = (183:366); 
a3 = 0.5*105;
b3 = 2*pi/365;
d3 = a3;
c3 = -93;
y3 = a3*sin(b3*(x3-c3))+d3;
%% Plot 
figure;
for i = 1:length(years)
    plot(1:length(rainfallData),rainfallData(:,i),'LineWidth',1.5)
    hold on
end
plot(x1,y1,'b--','LineWidth',1.5)
hold on
plot(x2,y2,'r--','LineWidth',1.5)
hold on
plot(x3,y3,'r--','LineWidth',1.5)
title('Comparison of Rainfall over 8 Years with Varying Approximations','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
xlim([0 366])
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')
lgd = legend({'2012','2013','2014','2015','2016','2017','2018','2019','2020','Cosine Approximation','Piecewise Approximation'},'Location','eastoutside');
lgd.FontSize = 20;
lgd.Interpreter = 'LaTeX';
%% Average over 8 years of data
rainfallTemp = rainfallData;
rainfallTemp(isnan(rainfallTemp))=0;
average = (rainfallTemp(:,1)+rainfallTemp(:,2)+rainfallTemp(:,3)+ ... 
    rainfallTemp(:,4)+rainfallTemp(:,5)+rainfallTemp(:,6)+rainfallTemp(:,7)...
    +rainfallTemp(:,8)+rainfallTemp(:,9))/9;

%%
figure;
plot(average)
title('Rainfall Averaged over 8 Years','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
xlim([0 366])
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')

max(average);

%% Fourier Transform
% t = 1:366; % time vector (days)
% rain = average; % rain data (mm)
% Ts = mean(diff(t));
% Fs = 1/Ts; 
% Fn = Fs/2;
% L = numel(t); 
% FTs = fft(rain - mean(rain))/L;
% Fv = linspace(0,1,fix(L/2)+1)*Fn;
% Iv = 1:numel(Fv);
% 
% figure 
% plot(t, abs(FTs(Iv))*2)
% grid 
% xlabel('')
% ylabel('')
%% 
% F = fft(average);
% plot(t,abs(F)/)
% xlim([0 366])
%% fourier formula
% a0 = 1/pi * int(f_data) * dx;
% an = 1/pi * int(f_data)*cos(n*x) * dx; 
% bn = 1/pi * int(f_data)*sin(n*x) * dx; 
% % 
% f = 0.5*a0+sum(an*cos(n*x))+sum(bn*sin(n*x));
% 


%% Averaging each month
% Month = [31 29 31 30 31 30 31 31 30 31 30 31];
% 
% January = zeros(Month(1),1);
% February = zeros(Month(2),1);
% March = zeros(Month(3),1);
% April = zeros(Month(4),1);
% May = zeros(Month(5),1);
% June = zeros(Month(6),1);
% July = zeros(Month(7),1);
% August = zeros(Month(8),1);
% September = zeros(Month(9),1);
% October = zeros(Month(10),1);
% November = zeros(Month(11),1);
% December = zeros(Month(12),1);
% 
% for i = 1:length(years)
%     January = January + rainfallTemp(1:Month(1),i);
%     February = February + rainfallTemp(sum(Month(1))+1:sum(Month(1:2)),i);
%     March = March + rainfallTemp(sum(Month(1:2)+1):sum(Month(1:3)),i);
%     April = April + rainfallTemp(sum(Month(1:3))+1:sum(Month(1:4)),i);
%     May = May + rainfallTemp(sum(Month(1:4))+1:sum(Month(1:5)),i);
%     June = June + rainfallTemp(sum(Month(1:5))+1:sum(Month(1:6)),i);
%     July = July + rainfallTemp(sum(Month(1:6))+1:sum(Month(1:7)),i);
%     August = August + rainfallTemp(sum(Month(1:7))+1:sum(Month(1:8)),i);
%     September = September + rainfallTemp(sum(Month(1:8))+1:sum(Month(1:9)),i);
%     October = October + rainfallTemp(sum(Month(1:9))+1:sum(Month(1:10)),i);
%     November = November + rainfallTemp(sum(Month(1:10))+1:sum(Month(1:11)),i);
%     December = December + rainfallTemp(sum(Month(1:11))+1:sum(Month(1:12)),i);
% end
% 
% Months = [January February March April May June July August September ...
%     October November December];

%% Fourier Calculation
t = linspace(1,366,366); 
N = 10;
%%
z = isnan(rainfallData(:,1));
index = find(z==1);
k = rainfallData(:,1)
k(z) = [];
t(z) = [];
%%
[a0, an, bn, s_approx, T] = trigFS(k', t, N);
%%
s_approx = [s_approx(1:index),NaN,s_approx(index+1:end)];
t = [t(1:index),NaN,t(index+1:end)];
figure
plot(t, s_approx)
below_0 = s_approx<0;
s_approx(below_0) =0;
figure
plot(t, s_approx)
%%
[a0, an, bn, s_approx, T] = trigFS(rainfallData(:,1)', t, N);
%%

for i = 1:size(rainfallData,2)
    [a0(i), an(i), bn(i), s_approx(i), T(i)] = trigFS(rainfallData(:,i)', t, N);
end

%% FUNCTION
function [a0, an, bn, s_approx, T] = trigFS(s_hinf, t, N)
    Ts = t(2) - t(1);
    T = t(end) - t(1) + Ts;
    f0 = 1/T;
    n = 1:N;
    a0 = 1/T*sum(s_hinf)*Ts; 
    an = 2/T *s_hinf * cos(2*pi*f0*t'*n) * Ts;
    bn = 2/T *s_hinf * sin(2*pi*f0*t'*n) * Ts;
    s_approx = a0 + an*cos(2*pi*f0*n'*t) + bn*sin(2*pi*f0*n'*t);
end






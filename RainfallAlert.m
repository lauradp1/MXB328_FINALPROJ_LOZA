%% Import Data
rainfall = readtable('MackayAlert.csv');
%% Remove unneccesary columns and convert to matrix
rainfall(:,[1,2,7,8]) = [];
rainfall = table2array(rainfall); %hello
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






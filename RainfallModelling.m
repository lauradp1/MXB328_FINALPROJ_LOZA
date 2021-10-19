%% Import Data
rainfall = readtable('MackayAlert.csv');
% Remove unneccesary columns and convert to matrix
rainfall(:,[1,2,7,8]) = [];
rainfall = table2array(rainfall);
rainfall2 = rainfall;
% Remove all other years outside of 2012-2020
rainfall = rainfall(rainfall(:,1)<2021 & rainfall(:,1)>=2011,:);
% Change format of data
years = unique(rainfall(:,1));
rainfallData = zeros(366,length(years));
for Average = 1:length(years)
    yearData = rainfall(rainfall(:,1)==years(Average),4);
    if length(yearData) ~= 366
        yearData = [yearData(1:59);NaN;yearData(60:end)];
    end
    rainfallData(:,Average) = yearData;
end
% Averaging the Data
rainfallTemp = rainfallData(:,[2:10]);
T = rainfallTemp';
average = mean(T,'omitnan')';
averageDaily = mean(average);
r_f = averageDaily*10^(-3);
% Replacing Nans with Averages
RainfallModel = rainfallData;
for j = 1:length(years)
    index = find(isnan(RainfallModel(:,j)));
    RainfallModel(index,j) = average(index);
end
%% Average Rainfall Plot
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
%% Cosine Approximation
rf = r_f;
syms t
%t = linspace(1,366,366);
cosineRain = @(t) rf + rf*cos(2*pi*t/366); % [mm/day] period of 366 days
figure
fplot(cosineRain,[0 366])

%% Fourier Calculation
% For average rainfall vector
t = linspace(1,366,366);
N = 10;
k = average;
[a0, an, bn, s_approx, T] = trigFS(k', t, N);
below_0 = s_approx<0;
s_approx(below_0) =0;
figure
plot(t, s_approx) %0.5, 0.6 
%% Average Plot with Fourier Approximation
t = linspace(1,366,366);
figure;
plot(t,average,'LineWidth',2)
hold on
plot(t,s_approx,'r--','LineWidth',2)
title('Rainfall Averaged from 2012 to 2020','FontSize',24,'Interpreter','LaTeX')
xlabel('Time (days)','FontSize',20,'Interpreter','LaTeX')
xlim([0 366])
ylabel('Rainfall (mm)','FontSize',20,'Interpreter','LaTeX')
lgd = legend('Averaged Data','Fourier Approximation');
lgd.Interpreter = 'LaTeX';
lgd.FontSize = 16;
%% Function 
%s = @(t) a0/2 + sum(an*cos(n*t) + bn*sin(n*t))
syms n 
s = @(t) a0/2 + symsym(an*cos(n*t) + bn*sin(n*t),1,inf);
%% 2012 Data
% For 2012 rainfall vector
t = linspace(1,366,366);
N = 50;
k = RainfallModel(:,2);
[a0, an, bn, s_approx, T] = trigFS(k', t, N);
below_0 = s_approx<0;
s_approx(below_0) =0;
figure
plot(t, s_approx)
%% Year of Flooding
% For 2011 rainfall vector
t = linspace(1,366,366);
N = 50;
k = RainfallModel(:,1);
[a0, an, bn, s_approx, T] = trigFS(k', t, N);
below_0 = s_approx<0;
s_approx(below_0) =0;
figure
plot(t, s_approx)
%% Flood -> Probability
Rain = RainfallModel(:,1) ~= 0;
rr = strfind(Rain', [1 1]);
dd = strfind(Rain', [0 0]);
T_flood = [length(dd)/365   1-length(dd)/365 ; 1-length(rr)/365 length(rr)/365];
% above matrix [dd, dr; rd, rr]
steadystate = T_flood^100; % probability of rain and no rain
%% Average -> Probability
RainA = average ~= 0;
rrA= strfind(RainA', [1 1]);
ddA = strfind(RainA', [0 0]);
T_floodA = [length(ddA)/365   1-length(ddA)/365 ; 1-length(rrA)/365 length(rrA)/365];
% above matrix [dd, dr; rd, rr]
steadystateA = T_floodA^100; % probability of rain and no rain
%%
% max(rainfallData(:,1))
% find(rainfallData(:,1) == max(rainfallData(:,1)))

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
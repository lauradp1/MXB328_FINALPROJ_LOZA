function [years, cosineRain, fourierRain, r_f, average, fourierFlood, cosineFlood, RainfallModel, steadystateFlood, steadystateAveraged, FourierAverage, FourierFlood] = RainfallModelling(rainfall)
% Climate Modelling Function 
% Takes an input of rainfall which is the data gathered from BOM
% Outputs:
% cosineRain: cosine approximation from averaged data
% fourierRain: fourier approximation from averaged data
% r_f: overall daily average (2012-2020)
% fourierFlood: fourier approximation from flood year data
% cosineFlood: cosine approximation from flood year data
% RainfallModel: overall matrix with complete and consistent data
% steadystateFlood: probability of rain for flood data
% steadystateAverage: probability of rain for averaged data
%% Remove unneccesary columns and convert to matrix
rainfall(:,[1,2,7,8]) = [];
rainfall = table2array(rainfall);
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
rainfallTemp = rainfallData(:,(2:10));
T = transpose(rainfallTemp);
average = mean(T,'omitnan')';
averageDaily = mean(average);
r_f = averageDaily;
% Replacing Nans with Averages
RainfallModel = rainfallData;
for j = 1:length(years)
    index = find(isnan(RainfallModel(:,j)));
    RainfallModel(index,j) = average(index);
end
%% Averaged Cosine Approximation
cosineRain = @(t) r_f + r_f*cos(2*pi*t/366); % [mm/day] period of 366 days
%% Average Fourier Approximation
t = linspace(1,366,366);
N = 50;
k = average;
[a0, an, bn, FourierAverage, T] = trigFS(k', t, N);
below_0 = FourierAverage < 0; 
FourierAverage(below_0) = 0;
fourierRain = @(t) a0 + sum(an(1:N) .*cos(2*pi*1/T.*(1:N).*t) + bn(1:N).*sin(2*pi*1/T.*(1:N).*t));
%% Flood Cosine Approximation
r_flood = mean(RainfallModel(:,1));
cosineFlood = @(t) r_flood + r_flood*cos(2*pi*t/366); % [mm/day] period of 366 days
%% Flood Fourier Approximation
t = linspace(1,366,366);
N = 50;
k = RainfallModel(:,1);
[a0_flood, an_flood, bn_flood, FourierFlood, T] = trigFS(k', t, N);
below_0 = FourierFlood < 0; 
FourierFlood(below_0) = 0;
fourierFlood = @(t) a0_flood + sum(an_flood(1:N) .*cos(2*pi*1/T.*(1:N).*t) + bn_flood(1:N).*sin(2*pi*1/T.*(1:N).*t));
%% Flood Data Probability
Rain = RainfallModel(:,1) ~= 0;
rr = strfind(Rain', [1 1]);
dd = strfind(Rain', [0 0]);
T_flood = [length(dd)/365   1-length(dd)/365 ; 1-length(rr)/365 length(rr)/365];
% above matrix [dd, dr; rd, rr]
steadystateFlood = T_flood^100; % probability of rain and no rain
%% Averaged Data Probability
RainA = average ~= 0;
rrA= strfind(RainA', [1 1]);
ddA = strfind(RainA', [0 0]);
T_floodA = [length(ddA)/365   1-length(ddA)/365 ; 1-length(rrA)/365 length(rrA)/365];
% above matrix [dd, dr; rd, rr]
steadystateAveraged = T_floodA^100; % probability of rain and no rain
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
function [cosineRain, fourierRain, r_f, average, fourierFlood, cosineFlood, RainfallModel] = RainfallModelling(rainfall)
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
r_f = averageDaily;
% Replacing Nans with Averages
RainfallModel = rainfallData;
for j = 1:length(years)
    index = find(isnan(RainfallModel(:,j)));
    RainfallModel(index,j) = average(index);
end
%% Cosine Approximation
rf = r_f;
cosineRain = @(t) rf + rf*cos(2*pi*t/366); % [mm/day] period of 366 days
%% Fourier Calculation
% For average rainfall vector
t = linspace(1,366,366);
N = 50;
k = average;
[a0, an, bn, sFourierAverage, T] = trigFS(k', t, N);
below_0 = sFourierAverage<0;
sFourierAverage(below_0)=0;
fourierRain = @(t) a0 + sum(an(1:N) .*cos(2*pi*1/T.*(1:N).*t) + bn(1:N).*sin(2*pi*1/T.*(1:N).*t));
% figure
% fplot(fun)
% xlim([0 366])
% ylim([0 14])
% Year of Flooding
%% For 2011 rainfall vector
r_flood = mean(RainfallModel(:,1));
cosineFlood = @(t) r_flood + r_flood*cos(2*pi*t/366); % [mm/day] period of 366 days
%% Fourier 
t = linspace(1,366,366);
N = 100;
k = RainfallModel(:,1);
[a0_flood, an_flood, bn_flood, s_flood, T] = trigFS(k', t, N);
fourierFlood = @(t) a0_flood + sum(an_flood(1:N) .*cos(2*pi*1/T.*(1:N).*t) + bn_flood(1:N).*sin(2*pi*1/T.*(1:N).*t));
% below_0 = s_flood<0;
% s_flood(below_0) =0;
% figure
% plot(t, s_flood)
% Flood -> Probability
% Rain = RainfallModel(:,1) ~= 0;
% rr = strfind(Rain', [1 1]);
% dd = strfind(Rain', [0 0]);
% T_flood = [length(dd)/365   1-length(dd)/365 ; 1-length(rr)/365 length(rr)/365];
% above matrix [dd, dr; rd, rr]
% steadystate = T_flood^100; % probability of rain and no rain
% Average -> Probability
% RainA = average ~= 0;
% rrA= strfind(RainA', [1 1]);
% ddA = strfind(RainA', [0 0]);
% T_floodA = [length(ddA)/365   1-length(ddA)/365 ; 1-length(rrA)/365 length(rrA)/365];
% above matrix [dd, dr; rd, rr]
%steadystateA = T_floodA^100; % probability of rain and no rain
% max(rainfallData(:,1))
% find(rainfallData(:,1) == max(rainfallData(:,1)))
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
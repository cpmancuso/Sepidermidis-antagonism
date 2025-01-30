function [fighandle,hExp, pExp] = plot_exponential_dist(data,fignum)

fighandle = figure(fignum);
clf(fignum)
fighandle.Position = [fighandle.Position(1) fighandle.Position(2) 500 200];

data = data';
% data = data(data>0)';

% Fit distributions
expFit = fitdist(data, 'Gamma');


% Generate the cumulative distribution functions (CDFs) 
[expCDF, xExp] = ecdf(data); % Empirical CDF for data
expTheoreticalCDF = cdf(expFit, xExp);


% Perform Kolmogorov-Smirnov tests
[hExp, pExp] = kstest(data, 'CDF', expFit);


% Plotting the fits
subplot(1,2,1)
histogram(data,[-1:2:60])
ylabel('Lineages')
xlabel('Number of Other Lineages Antagonized')

subplot(1,2,2)
plot(xExp, expCDF, 'k-', 'Marker','.')
hold on
plot(xExp, expTheoreticalCDF, 'r--');
title(['Exponential, KS p = ' num2str(pExp,'%.3f')]);
ylabel('Density')
xlabel('Number of Other Lineages Antagonized')

set(gca,'box','off')
pbaspect([1 1 1])
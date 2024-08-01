%% Dynamical Brain Mapping Group - DYNAMAP
%  Institut de Neurosciences des Systemes-  INS
%  Faculte de Medecine, Aix-Marseille Universite
%  27, Boulevard Jean Moulin, 13005 Marseille, France
%  https://ins-amu.fr/dynamap

%Authors: Damian Dellavale (dellavaledamian@gmail.com
%                           hector.dellavale-clara@univ-amu.fr, 
%                           dellavale@cab.cnea.gov.ar).

%Project: Local false discovery rate analysis.

%Date: 10/06/2021

%% Description.

%In this function, we implement the Local False Discovery Rate (LFDR) method,
%which is a Bayesian approach to control type I errors (false positives)
%in defining a significance threshold corrected for (many) multiple comparisons.

%NOTE: The LFDR method does not necessarily have to be applied on P values or z values:
%"It is convenient, although not necessary, to work with z-values instead of the Yi’s or Pi’s"
%Ref:
%/Papers/statistical analysis/FalseDiscoveryRate/Efron2004.pdf

%% Tree of dependencies.

%function_LFDR_v1.m
%

%% Changes from previous versions.

%None.

%% How to test this script.

%\matlab_functions\statistics\lfdr\_testing\example_lfdr_v1.m

%% References.

%Bibliography used to implement this function:
%
%/Papers/statistical analysis/FalseDiscoveryRate
%Efron2004.pdf
%Efron_Size, power and false discovery rates.pdf
%Efron_2005LocalFDR.pdf
%
%https://en.wikipedia.org/wiki/Bayes'_theorem
%
%https://en.wikipedia.org/wiki/Probability_distribution
%https://en.wikipedia.org/wiki/Probability_density_function
%https://en.wikipedia.org/wiki/Cumulative_distribution_function
%https://en.wikipedia.org/wiki/Kernel_(statistics)
%
%Half-normal_distribution is the probability distribution of y = |x|, 
%where x is a random variable following an ordinary normal distribution.
%https://en.wikipedia.org/wiki/Folded_normal_distribution
%https://en.wikipedia.org/wiki/Half-normal_distribution
%https://fr.mathworks.com/help/stats/prob.halfnormaldistribution.html
%
%https://fr.mathworks.com/help/stats/ksdensity.html
%https://fr.mathworks.com/help/stats/fitdist.html
%https://fr.mathworks.com/help/stats/prob.normaldistribution.pdf.html
%https://fr.mathworks.com/help/stats/makedist.html
%https://fr.mathworks.com/help/stats/normpdf.html
%https://fr.mathworks.com/help/stats/fitgmdist.html
%https://fr.mathworks.com/help/stats/normfit.html
%
%https://en.wikipedia.org/wiki/Interquartile_range
%https://fr.mathworks.com/help/stats/prob.normaldistribution.iqr.html
%https://fr.mathworks.com/help/stats/quantile.html

%Scripts of reference used to implement this function:
%
%https://gitlab-dynamap.timone.univ-amu.fr/epitools/dyn_toolbox/-/blob/feature/lfdr/stats/dyn_lfdr.m
%/delphos4bids/fast_lfdr.m
%/DELPHOS-Delphos1.5/Delphos.m
%
%/projects/neuro/PDbiomarkers/animalModels/IFIBIO_UBA/IFIBIO_14d12d2018/damian/ged_bandPassFiltering/concatSamples/
%ged_twoGroupAnalysis_individualGED_v1.m
%ged_twoGroupAnalysis_individualGED_v2.m

%Official LFDR package in R:
%
%http://www.strimmerlab.org/notes/fdr.html
%
%fdrtool: a versatile R package for estimating local and tail area-based false discovery rates
%https://academic.oup.com/bioinformatics/article/24/12/1461/196272
%
%https://cran.r-project.org/web/packages/locfdr/
%https://cran.r-project.org/web/packages/locfdr/locfdr.pdf
%https://cran.r-project.org/web/packages/locfdr/vignettes/locfdr-example.pdf
%https://rdrr.io/cran/locfdr/man/locfdr.html

%% Main function.

function [LFDRout] = function_LFDR_v1(X, varargin)
%==========================================================================
%Inputs:
%
%X  -> Input signal (array: 1 x samples = 1 x Nsamples).
%
%LFDRcfg -> Configuration structure with fields:
%           'histBinMethod'         -> 
%           'histNormalization'     -> 
%
%           'mixtureDistFunction'   -> 
%           'mixtureDistName'       -> 
%           'mixtureDistParam'      -> 
%
%           'theoreticalDistName'   -> 
%           'theoreticalDistParam'  -> 
%
%           'empiricalDistRange'    -> 
%           'empiricalDistFunction' -> 
%           'empiricalDistName'     -> 
%           'empiricalDistParam'    -> 
%
%           'lfdr_threshold'        -> 
%
%           'plotFlag'              -> 

%Outputs:
%
%LFDRout -> Output structure with fields:
%
%           'theoThreshold_low'  -> 
%           'theoThreshold_high' -> 
%
%           'empThreshold_low'   -> 
%           'empThreshold_high'  -> 
%
%           'LFDRcfg'            -> 
%==========================================================================

%Check the input arguments ------------------------------------------------

args = checkInputs(X, varargin);

%--------------------------------------------------------------------------

%Default values of the outputs --------------------------------------------
%--------------------------------------------------------------------------

%Parameters ---------------------------------------------------------------
dimRows = 1;
dimCols = 2;

dimSamples = 2; %columns => samples.

Nsamples = size(X,dimSamples); %Number of samples.
%--------------------------------------------------------------------------

%Compute the Local False Discovery Rate (LFDR) method ---------------------

%---
%Compute the histogram of the X values.

%Example of the input parameters configuration:
%LFDRcfg.histBinMethod = 'auto';
%LFDRcfg.histNormalization = 'pdf';

[counts, edges] = histcounts(X,...
                  'BinMethod', args.histBinMethod,...
                  'Normalization', args.histNormalization);

% %Equivalent to 'pdf' normalization.
% [counts, edges] = histcounts(X, 'histBinMethod', args.histBinMethod);              
% counts = counts ./ (Nsamples * (edges(2:end)-edges(1:end-1)));              

%Compute the central value of the bins.
Xbins = mean([edges(1:end-1); edges(2:end)], dimRows);
%---

%---
%Compute the MIXTURE probability density function of X.

switch lower(args.mixtureDistFunction)

    case 'ksdensity'
        
        %Example of the input parameters configuration:
        %LFDRcfg.mixtureDistFunction = 'ksdensity';
        %LFDRcfg.mixtureDistName = '';
        %LFDRcfg.mixtureDistParam = {'Function', 'pdf'};        
        
        %Kernel smoothing function estimate.
        [f, Xbins] = ksdensity(X, Xbins, args.mixtureDistParam{:});

    case 'fitdist'
        
        %Example of the input parameters configuration:
        %LFDRcfg.mixtureDistFunction = 'fitdist';
        %LFDRcfg.mixtureDistName = 'Kernel';
        %LFDRcfg.mixtureDistParam = {'Kernel', 'epanechnikov'};        
        
        %Fit a probability distribution to the X values.
        pdObject = fitdist(X.', args.mixtureDistName, args.mixtureDistParam{:});
        f = pdf(pdObject, Xbins);

    otherwise
        
        error('MATLAB:localFDR',['Mixture distribution function not identified.']);

end %switch lower(args.mixtureDistFunction)
%---

%---
%Compute the probability density function of X under the THEORETICAL null hypothesis.

switch lower(args.theoreticalDistFunction)
    
    case 'makedist'
        
        %Example of the input parameters configuration:
        %LFDRcfg.theoreticalDistFunction = 'makedist';
        %LFDRcfg.theoreticalDistName = 'Normal';
        %LFDRcfg.theoreticalDistParam = {'mu', 0, 'sigma', 1};
        
        %Custom probability density function.
        pdObject = makedist(args.theoreticalDistName, args.theoreticalDistParam{:});
        f0_theo = pdf(pdObject, Xbins);

        % %Another way to compute a normal probability density function.
        % f0_theo = normpdf(Xbins, args.theoreticalDistParam{:});

    case 'fitdist'

        %Example of the input parameters configuration:
        %LFDRcfg.theoreticalDistFunction = 'fitdist';
        %LFDRcfg.theoreticalDistName = 'Normal';
        %LFDRcfg.theoreticalDistParam = {};

        %Fit a probability distribution to the X values.
        pdObject = fitdist(X.', args.theoreticalDistName, args.theoreticalDistParam{:});
        f0_theo = pdf(pdObject, Xbins);

    case 'normfit'

        %Example of the input parameters configuration:
        %LFDRcfg.theoreticalDistFunction = 'normfit';
        %LFDRcfg.theoreticalDistName = 'Normal';
        %LFDRcfg.theoreticalDistParam = {};
            
        %Fit a normal probability distribution to the X values.
        %Normal parameter estimates.
        [mu, sigma] = normfit(X);
        args.theoreticalDistParam = {'mu', mu, 'sigma', sigma};
        
        %Normal probability density function.
        pdObject = makedist(args.theoreticalDistName, args.theoreticalDistParam{:});
        f0_theo = pdf(pdObject, Xbins);
        %
        %f0_theo = normpdf(Xbins, args.theoreticalDistParam{:});

    otherwise
        
        error('MATLAB:localFDR',['Theoretical distribution function not identified.']);

end %switch lower(args.theoreticalDistFunction)
%---

%---
%Compute the probability density function of X under the EMPIRICAL null hypothesis.
%IMPORTANT: Situations occur in which going from the theoretical to empirical
%null hypothesis totally change the findings of significance. 

%Compute the intervals to fit the probability density function of X under the EMPIRICAL null hypothesis.
if strcmpi(args.empiricalDistRangeMethod, 'auto')
    %Explore a range of intervals to fit the EMPIRICAL probability density function of X.
    empiricalDistFittingInt =...
    [linspace(args.empiricalDistRange(1), mean(X), args.empiricalDistNfitInt).',...
     flipud(linspace(mean(X), args.empiricalDistRange(end), args.empiricalDistNfitInt).')];
else
    %Use a single interval to fit the EMPIRICAL probability density function of X.
    empiricalDistFittingInt = args.empiricalDistRange;
end %if strcmpi(args.empiricalDistRangeMethod, 'auto')

%Memory pre-allocation to speed up the loop.         
f0_emp_ = NaN(size(empiricalDistFittingInt,dimRows),length(f));
Erms = NaN(1,size(empiricalDistFittingInt,dimRows));  
for ii=1:+1:size(empiricalDistFittingInt,1) %Loop over the fitting intervals.

%Compute the Xbins_ bins inside the current fitting interval.  
ind_Xbins = Xbins > empiricalDistFittingInt(ii,1) & Xbins < empiricalDistFittingInt(ii,2);
Xbins_ = Xbins(ind_Xbins);

if ~isempty(Xbins_)

%Compute the signal values X_ for the bins Xbins_ corresponding to the current fitting interval.   
ind_X = X > Xbins_(1) & X < Xbins_(end);
X_ = X(ind_X);

%Fit the probability density function for the current fitting interval.
switch lower(args.empiricalDistFunction)
    
    case 'fitdist'
        
        %Example of the input parameters configuration:
        %LFDRcfg.empiricalDistFunction = 'fitdist';
        %LFDRcfg.empiricalDistName = 'Normal';
        %LFDRcfg.empiricalDistParam = {};
        
        %Fit a probability distribution to the X_ values.
        try
            pdObject = fitdist(X_.', args.empiricalDistName, args.empiricalDistParam{:});
            f0_emp_(ii,:) = pdf(pdObject, Xbins);
        catch
            continue,
        end %try
            
    case 'normfit'
        
        %Example of the input parameters configuration:
        %LFDRcfg.empiricalDistFunction = 'normfit';
        %LFDRcfg.empiricalDistName = 'Normal';
        %LFDRcfg.empiricalDistParam = {};        
        
        %Fit a normal probability distribution to the X_ values.
        try
            %Normal parameter estimates.
            [mu_, sigma_] = normfit(X_);
            args.empiricalDistParam = {'mu', mu_, 'sigma', sigma_};
        
            %Normal probability density function.
            pdObject = makedist(args.empiricalDistName, args.empiricalDistParam{:});
            f0_emp_(ii,:) = pdf(pdObject, Xbins);
            %
            %f0_emp_(ii,:) = normpdf(Xbins, args.empiricalDistParam{:});
        catch
            continue,
        end %try
        
    otherwise
        
        error('MATLAB:localFDR',['Empirical distribution function not identified.']);

end %switch lower(args.empiricalDistFunction)

%Define the interval to compute the RMS error between: 
%the MIXTURE probability density function (f) and
%the probability density function under the EMPIRICAL null hypothesis (f0_emp_).
%indErms = true(size(f)); %The whole interval of bins.
indErms = ind_Xbins; %The bins inside the current fitting interval.

%Compute the RMS error between: 
%the MIXTURE probability density function (f) and
%the probability density function under the EMPIRICAL null hypothesis (f0_emp_).
%
%RMS error using the non-weighted difference.
%Erms(ii) = sqrt( sum((f(indErms) - f0_emp_(ii,indErms)).^2)/length(f(indErms)) );
%
%RMS error using the weighted difference to assign less importance to the 
%difference in the tails of the distributions.
Erms(ii) = sqrt( sum(f(indErms).*(f(indErms) - f0_emp_(ii,indErms)).^2)/length(f(indErms)) );

end %if ~isempty(Xbins_)

end %Loop over the fitting intervals.

%Compute the coordinates for the minimum RMS error. 
[valMin, indMin] = min(Erms);

%Compute the Xbins_ bins corresponding to the fitting interval 
%producing the minimum RMS error.
ind_Xbins = Xbins > empiricalDistFittingInt(indMin,1) & Xbins < empiricalDistFittingInt(indMin,2);
Xbins_ = Xbins(ind_Xbins);

%Compute the fitting interval corresponding to the fitting interval 
%producing the minimum RMS error.
args.empiricalDistFittingInt = empiricalDistFittingInt(indMin,:);

%Compute the the probability density function under the EMPIRICAL null
%hypothesis (f0_emp), corresponding to the fitting interval producing the
%minimum RMS error.
f0_emp = f0_emp_(indMin,:);
%---

%---
%Compute the false discovery rate using the Bayes theorem: 
%a posteriori probability of being in the uninteresting class given a X value.

lfdr_theo = f0_theo ./ f;
lfdr_emp = f0_emp ./ f;
%---

%---
%Compute the thresholds for X values corresponding to a maximum false 
%discovery rate value (i.e. maximum rate of false positive X values).

%Thresholds for the THEORETICAL PDF.
if sum((lfdr_theo < args.lfdr_threshold) & (Xbins < 0))
    theoThreshold_low = Xbins( find((lfdr_theo < args.lfdr_threshold) & (Xbins < 0), 1, 'last') );
else
    theoThreshold_low = Xbins(1);
end %if sum((lfdr_theo < args.lfdr_threshold) & (Xbins < 0))

if sum((lfdr_theo < args.lfdr_threshold) & (Xbins > 0))
    theoThreshold_high = Xbins( find((lfdr_theo < args.lfdr_threshold) & (Xbins > 0), 1, 'first') );
else
    theoThreshold_high = Xbins(end);
end %if sum((lfdr_theo < args.lfdr_threshold) & (Xbins > 0))

%Thresholds for the EMPIRICAL PDF.
if sum((lfdr_emp < args.lfdr_threshold) & (Xbins < 0))
    empThreshold_low = Xbins( find((lfdr_emp < args.lfdr_threshold) & (Xbins < 0), 1, 'last') );
else
    empThreshold_low = Xbins(1);
end %if sum((lfdr_emp < args.lfdr_threshold) & (Xbins < 0))

if sum((lfdr_emp < args.lfdr_threshold) & (Xbins > 0))
    empThreshold_high = Xbins( find((lfdr_emp < args.lfdr_threshold) & (Xbins > 0), 1, 'first') );
else
    empThreshold_high = Xbins(end);
end %if sum((lfdr_emp < args.lfdr_threshold) & (Xbins > 0))
%---

if args.plotFlag %Plot the results.

%---    
%Compute the fit intervals for valid RMS error values, i.e. intervals for
%which Xbins_ is non empty.
empiricalDistFittingInt = empiricalDistFittingInt(1:length(Erms),:);

%SUBPLOT 1: Show the probability density functions.
figure, 
ax(1) = subplot(3,1,1);
hold on,
ph(1) = plot(Xbins,f,'-b','LineWidth', 2);
ph(2) = plot(Xbins,f0_emp_(1,:),'--','LineWidth',1);
plot(Xbins,f0_emp_,'--','LineWidth',1);

box on, grid off, axis tight, haxes = gca;

set(gca,'xscale','linear');
%set(gca,'yscale','linear');
%set(gca,'xscale','log');
set(gca,'yscale','log');

%set(gca,'XTickLabel',{'';''});

xlabel('X [arb. units]','interpreter','none')
ylabel(['Probability density function (PDF)'],'interpreter','none')

%Legends.
legend(ph,...
       'Mixture PDF',...
       'Empirical PDF',...
       'location','best','interpreter','none');
%legend('boxoff')
legend('boxon') 

%Clear the plot handle.
clear ph

%SUBPLOT 2: Show the resulting Empirical PDF.
ax(2) = subplot(3,1,2);
hold on,
ph(1) = plot(Xbins,f,'-b','LineWidth',2);
ph(2) = plot(Xbins,f0_emp_(indMin,:),'-r','LineWidth',2);

box on, grid off, axis tight, haxes = gca;

% for ii=1:+1:size(empiricalDistFittingInt,1) %Loop over the fitting intervals.
% plot([empiricalDistFittingInt(ii,1), empiricalDistFittingInt(ii,1)],...
%      [haxes.YLim(1), haxes.YLim(end)],...
%      'Color',[ii/size(empiricalDistFittingInt,1), ii/size(empiricalDistFittingInt,1), ii/size(empiricalDistFittingInt,1)],...
%      'LineWidth',1);
% plot([empiricalDistFittingInt(ii,2), empiricalDistFittingInt(ii,2)],...
%      [haxes.YLim(1), haxes.YLim(end)],...
%      'Color',[ii/size(empiricalDistFittingInt,1), ii/size(empiricalDistFittingInt,1), ii/size(empiricalDistFittingInt,1)],...
%      'LineWidth',1)
% end %Loop over the fitting intervals.

ii = 1;
ph(3) = plot([empiricalDistFittingInt(ii,1), empiricalDistFittingInt(ii,1)],[haxes.YLim(1), haxes.YLim(end)], '--k', 'LineWidth', 1);
ph(3) = plot([empiricalDistFittingInt(ii,2), empiricalDistFittingInt(ii,2)],[haxes.YLim(1), haxes.YLim(end)], '--k', 'LineWidth', 1);
% ii = size(empiricalDistFittingInt,1);
% ph(3) = plot([empiricalDistFittingInt(ii,1), empiricalDistFittingInt(ii,1)],[haxes.YLim(1), haxes.YLim(end)],'--k','LineWidth',1);
% ph(3) = plot([empiricalDistFittingInt(ii,2), empiricalDistFittingInt(ii,2)],[haxes.YLim(1), haxes.YLim(end)],'--k','LineWidth',1);

ii = indMin;
ph(4) = plot([empiricalDistFittingInt(ii,1), empiricalDistFittingInt(ii,1)],[haxes.YLim(1), haxes.YLim(end)],'--r','LineWidth',1);
ph(4) = plot([empiricalDistFittingInt(ii,2), empiricalDistFittingInt(ii,2)],[haxes.YLim(1), haxes.YLim(end)],'--r','LineWidth',1);

%set(gca,'XTickLabel',{'';''});

set(gca,'xscale','linear');
set(gca,'yscale','linear');
%set(gca,'xscale','log');
%set(gca,'yscale','log');

xlabel('X [arb. units]','interpreter','none')
ylabel(['Probability density function (PDF)'],'interpreter','none')

%SUBPLOT 3: Show the RMS error.
subplot(3,1,3);
hold on,
plot(empiricalDistFittingInt(:,2)-empiricalDistFittingInt(:,1),Erms,':ob')
plot(empiricalDistFittingInt(indMin,2)-empiricalDistFittingInt(indMin,1),Erms(indMin),'+r','MarkerSize',20)

box on, grid off, axis tight, haxes = gca;

%set(gca,'XTickLabel',{'';''});

set(gca,'xscale','linear');
%set(gca,'yscale','linear');
%set(gca,'xscale','log');
set(gca,'yscale','log');

xlabel('Fitting interval around the mean of X [arb. units]','interpreter','none')
%ylabel(['RMS error between the Mixture PDF and the fitted Empirical PDF'],'interpreter','none')
ylabel(['RMS error'],'interpreter','none')

%Link the scale for the x-axis between the plots.
linkaxes(ax, 'x');

%Legends.
legend(ph,...
       'Mixture PDF',...
       'Empirical PDF',...
       'Maximum allowed fitting interval',...
       'Fitting interval producing the minimum RMS error',...
       'location','best','interpreter','none');
%legend('boxoff')
legend('boxon') 

%Clear the plot and axis handles.
clear ph ax
%---    
    
%---    
%SUBPLOT 1: Show the probability density functions.
%rgbBLUE = [0,0,1];
%rgbDARKBLUE  = [0,0,0.7];
%rgbLIGHTBLUE  = [0.3,0.3,1];
rgbVERYLIGHTBLUE  = [0.5,0.5,1];
rgbWHITE = [1,1,1];

figure, 
ax(1) = subplot(2,1,1);
hold on,
ph(1) = bar(Xbins, counts, 'FaceColor', rgbVERYLIGHTBLUE, 'EdgeColor', rgbWHITE);
ph(2) = plot(Xbins, f, '-r', 'LineWidth', 2);
ph(3) = plot(Xbins, f0_theo, ':k', 'LineWidth', 2);
ph(4) = plot(Xbins, f0_emp, '--b', 'LineWidth', 2);
%ph(1) = bar(Xbins, log10(counts), 'FaceColor', rgbVERYLIGHTBLUE, 'EdgeColor', rgbWHITE);
%ph(2) = plot(Xbins, log10(f), '-r', 'LineWidth', 2);
%ph(3) = plot(Xbins, log10(f0_theo), ':k', 'LineWidth', 2);
%ph(4) = plot(Xbins, log10(f0_emp), '--b', 'LineWidth', 2);

box on, grid off, axis tight, haxes = gca;

%Vertical markers for the range of X values used to fit the empirical probability density function.
ph(5) = plot([Xbins_(1), Xbins_(1)], [haxes.YLim(1), haxes.YLim(2)], '--b', 'LineWidth', 1);
ph(5) = plot([Xbins_(end), Xbins_(end)], [haxes.YLim(1), haxes.YLim(2)], '--b', 'LineWidth', 1);
%ph(5) = plot([Xbins_(1), Xbins_(1)], log10([haxes.YLim(1), haxes.YLim(2)]), '--b', 'LineWidth', 1);
%ph(5) = plot([Xbins_(end), Xbins_(end)], log10([haxes.YLim(1), haxes.YLim(2)]), '--b', 'LineWidth', 1);

%set(gca,'XTickLabel',{'';''});

xlabel('X [arb. units]','interpreter','none')
ylabel(['Probability density function (PDF)'],'interpreter','none')

%Legends.
legend(ph,...
       'Histogram',...
       'Mixture PDF',...
       'Theoretical PDF',...
       'Empirical PDF',...
       'Interval used to fit the empirical PDF',...
       'location','best','interpreter','none');
%legend('boxoff')
legend('boxon') 

%Clear the plot handle.
clear ph

%SUBPLOT 2: Show the local false discovery rate.
ax(2) = subplot(2,1,2);
hold on,
ph(1) = plot(Xbins, log10(lfdr_theo), ':k', 'LineWidth', 2);
ph(2) = plot(Xbins, log10(lfdr_emp), '--b', 'LineWidth', 2);

%Plot the data points.
ph(3) = plot(X, log10(args.lfdr_threshold * ones(Nsamples,1)), 'x', 'MarkerSize', 10, 'Color', rgbVERYLIGHTBLUE);

box on, grid off, axis tight, haxes = gca;

%Horizontal marker showing the local FDR threshold.
ph(4) = plot([haxes.XLim(1), haxes.XLim(2)], log10([args.lfdr_threshold, args.lfdr_threshold]), '-k', 'LineWidth', 1);

%Vertical markers for the range of X values corresponding to the theoretical local FDR.
ph(5) = plot([theoThreshold_low, theoThreshold_low], [haxes.YLim(1), haxes.YLim(2)], ':k', 'LineWidth', 1);
ph(5) = plot([theoThreshold_high, theoThreshold_high], [haxes.YLim(1), haxes.YLim(2)], ':k', 'LineWidth', 1);

%Vertical markers for the range of X values corresponding to the empirical local FDR.
ph(6) = plot([empThreshold_low, empThreshold_low], [haxes.YLim(1), haxes.YLim(2)], '--b', 'LineWidth', 1);
ph(6) = plot([empThreshold_high, empThreshold_high], [haxes.YLim(1), haxes.YLim(2)], '--b', 'LineWidth', 1);

xlabel('X [arb. units]','interpreter','none')
ylabel(['log local FDR'],'interpreter','none')

%Link the scale for the x-axis between the plots.
linkaxes(ax, 'x');

%Legends.
legend(ph,...
       'Theoretical local FDR',...
       'Empirical local FDR',...
       'Data points',...
       'local FDR threshold',...
       'Threshold crossing for the theoretical PDF',...
       'Threshold crossing for the empirical PDF',...
       'location','best','interpreter','none');
%legend('boxoff')
legend('boxon') 

%Clear the plot and axis handles.
clear ph ax
%---

end %if args.plotFlag %Plot the results.

%--------------------------------------------------------------------------

%Synthesize the output structure ------------------------------------------

LFDRout = struct(...
'lfdr_emp', lfdr_emp,...
'empThreshold_low', empThreshold_low,...
'empThreshold_high', empThreshold_high,...
'lfdr_theo', lfdr_theo,...
'theoThreshold_low', theoThreshold_low,...
'theoThreshold_high', theoThreshold_high,...
'LFDRcfg', args);

%--------------------------------------------------------------------------

end %Main function.

%% checkInputs function.

function args = checkInputs(X, varargin)

%---
%Create the input parser object.
p = inputParser();

%Set the name of function for error message.
p.FunctionName = mfilename; %Use the name of the current running script.

%Set case-sensitive matches.
p.CaseSensitive = true;

%Throws an error if an input argument name does not match one defined in the input parser scheme. 
p.KeepUnmatched = true;

%Set Expanding structures into separate inputs, where each field name
%corresponds to an input parameter name. 
p.StructExpand = true;
%---

%---
%Define general constrains on the input arguments.

%validScalar = @(x) (isnumeric(x) && isscalar(x) && ~isnan(x) && ~isempty(x));
validScalar = @(x) isnumeric(x) && isscalar(x) && ~isnan(x);
validArray = @(x) isnumeric(x) && ~any(isnan(x(:)));
validString = @(x) isstring(x) || ischar(x);

%Add required, positional argument into input parser scheme.
validX = validArray;
p.addRequired('X', validX);

%Add optional name-value pair argument into input parser scheme.
%
%HISTOGRAM.
%
validHistBinMethod = validString;
%Default option 'auto': A Matlab algorithm chooses a bin width to cover the
%data range and reveal the shape of the underlying distribution.
defaultHistBinMethod = 'auto';
p.addParameter('histBinMethod', defaultHistBinMethod, validHistBinMethod);
%
validHistNormalization = validString;
defaultHistNormalization = 'pdf'; %Probability density function estimate.
p.addParameter('histNormalization', defaultHistNormalization, validHistNormalization);
%
%MIXTURE probability density function.
%
validMixtureDistFunction = validString;
defaultMixtureDistFunction = 'ksdensity'; %Kernel smoothing function estimate.
p.addParameter('mixtureDistFunction', defaultMixtureDistFunction, validMixtureDistFunction);
%
validMixtureDistName = validString;
%No distribution name is required to compute the default Kernel smoothing function estimate.
defaultMixtureDistName = '';
p.addParameter('mixtureDistName', defaultMixtureDistName, validMixtureDistName);
%
validMixtureDistParam = @(x) iscell(x); %@(x) (iscell(x) && ~isempty(x));
defaultMixtureDistParam = {'Function', 'pdf'}; %Probability density function.
p.addParameter('mixtureDistParam', defaultMixtureDistParam, validMixtureDistParam);
%
%THEORETICAL probability density function.
%
validTheoreticalDistFunction = validString;
defaultTheoreticalDistFunction = 'makedist'; %Create probability distribution object.
p.addParameter('theoreticalDistFunction', defaultTheoreticalDistFunction, validTheoreticalDistFunction);
%
validTheoreticalDistName = validString;
defaultTheoreticalDistName = 'Normal'; %Normal distribution.
p.addParameter('theoreticalDistName', defaultTheoreticalDistName, validTheoreticalDistName);
%
validTheoreticalDistParam = @(x) iscell(x); %@(x) (iscell(x) && ~isempty(x));
defaultTheoreticalDistParam = {'mu', 0, 'sigma', 1}; %Standard normal distribution.
p.addParameter('theoreticalDistParam', defaultTheoreticalDistParam, validTheoreticalDistParam);
%
%EMPIRICAL probability density function.
%
validEmpiricalDistRange = validArray;
%defaultEmpiricalDistRange = [quantile(X,1/4), quantile(X,3/4)]; %Interquartile range.
%defaultEmpiricalDistRange = [median(X) - iqr(X)/2, median(X) + iqr(X)/2]; %Interquartile range.
%defaultEmpiricalDistRange = [quantile(X,1/4) - 1.5*iqr(X), 1.5*iqr(X) + quantile(X,3/4)];
%defaultEmpiricalDistRange = [quantile(X,1/4) - 3*iqr(X), 3*iqr(X) + quantile(X,3/4)];
defaultEmpiricalDistRange = [quantile(X,1/4) - 5*iqr(X), 5*iqr(X) + quantile(X,3/4)];
%defaultEmpiricalDistRange = Inf * [-1, 1];
p.addParameter('empiricalDistRange', defaultEmpiricalDistRange, validEmpiricalDistRange);
%Refs:
%https://en.wikipedia.org/wiki/Interquartile_range
%
validEmpiricalDistRangeMethod = validString;
defaultEmpiricalDistRangeMethod = '';
p.addParameter('empiricalDistRangeMethod', defaultEmpiricalDistRangeMethod, validEmpiricalDistRangeMethod);
%
validEmpiricalDistNfitInt = validScalar;
defaultEmpiricalDistNfitInt = 100;
p.addParameter('empiricalDistNfitInt', defaultEmpiricalDistNfitInt, validEmpiricalDistNfitInt);
%
validEmpiricalDistFunction = validString;
defaultEmpiricalDistFunction = 'fitdist'; %Fit probability distribution object to data.
p.addParameter('empiricalDistFunction', defaultEmpiricalDistFunction, validEmpiricalDistFunction);
%
validEmpiricalDistName = validString;
defaultEmpiricalDistName = 'Normal'; %Normal distribution.
p.addParameter('empiricalDistName', defaultEmpiricalDistName, validEmpiricalDistName);
%
validEmpiricalDistParam = @(x) iscell(x); %@(x) (iscell(x) && ~isempty(x));
defaultEmpiricalDistParam = {}; %The parameters for the normal distribution are computed inside the function.
p.addParameter('empiricalDistParam', defaultEmpiricalDistParam, validEmpiricalDistParam);
%
%Threshold for the local FDR.
%
validLFDR_threshold = validScalar;
defaultLFDR_threshold = 0.1; % 10% of false positives.
p.addParameter('lfdr_threshold', defaultLFDR_threshold, validLFDR_threshold);
%
%Plot parameters.
%
validPlotFlag = @(x) islogical(x);
defaultPlotFlag = false; %No plot.
p.addParameter('plotFlag', defaultPlotFlag, validPlotFlag);
%
%NOTE: We desided not to use the positional argument format.
%For instance:
%
% %Add optional, positional argument into input parser scheme.
% validFs = validScalar;
% defaultFs = 1000; %[Hz] Sampling rate.
% p.addOptional('fs', defaultFs, validFs);

%ADD MORE INPUT PARAMETERS HERE...
%---

%---
%Parse the input arguments.
%This code line supports the input arguments be included in a single 
%structure (i.e. DDcfg), or in the name-value pair format.
p.parse(X, varargin{1,1}{:});

%Return the valid arguments.
args = p.Results;
%---

%---
%Define specific constrains on the input arguments.

%Check the input time series is a row array.
assert(size(X,2)>size(X,1),...
       'The input time series must be a row array.'); 

%If the range of X values to compute the empirical distribution is empty, assign to it the default value.
if isempty(args.empiricalDistRange)
    %args.empiricalDistRange = [quantile(X,1/4), quantile(X,3/4)]; %Interquartile range.
    %args.empiricalDistRange = [median(X) - iqr(X)/2, median(X) + iqr(X)/2]; %Interquartile range.
    %args.empiricalDistRange = [quantile(X,1/4) - 1.5*iqr(X), 1.5*iqr(X) + quantile(X,3/4)];
    %args.empiricalDistRange = [quantile(X,1/4) - 3*iqr(X), 3*iqr(X) + quantile(X,3/4)];
    args.empiricalDistRange = [quantile(X,1/4) - 5*iqr(X), 5*iqr(X) + quantile(X,3/4)];
    %args.empiricalDistRange = Inf * [-1, 1];
end %if isempty(args.empiricalDistRange)

%ADD MORE SPECIFIC CONSTRAINS HERE...
%---

%Refs:
%https://www.mathworks.com/help/matlab/ref/nargin.html
%https://www.mathworks.com/help/matlab/ref/varargin.html
%https://fr.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html
%https://www.mathworks.com/help/matlab/ref/inputparser.html            
            
end %checkInputs function.

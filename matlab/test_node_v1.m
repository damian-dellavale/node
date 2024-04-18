%%

%Consejo Nacional de Investigaciones Científicas y Técnicas (CONICET)
%Departamento de Física Médica (DFM)
%Instituto de Nanociencia y Nanotecnología (INN)
%Centro Atómico Bariloche (CAB)
%R8402AGP, San Carlos de Bariloche, Río Negro, Argentina
%https://inn-cnea.conicet.gov.ar/

% Dynamical Brain Mapping Group (DYNAMAP)
% Institut de Neurosciences des Systemes (INS)
% Faculte de Medecine, Aix-Marseille Universite
% 27, Boulevard Jean Moulin, 13005 Marseille, France
% https://ins-amu.fr/dynamap

% Authors: Damian Dellavale (dellavale@cab.cnea.gov.ar,
%                            dellavaledamian@gmail.com,
%                            https://github.com/damian-dellavale.github.io).

% Project: Nested Outlier Detection (NODE) algorithm.

% Date: 04/04/2023

%% Description.

%This is a test script for detecting events using the NODE algorithm (function_node_v1.m).
%For this, we use human intracranial EEG (iEEG) data. The example iEEG dataset
%was acquired at the Medical Center of the University of California, Irvine.
%The Office for the Protection of Human Subjects of the University of 
%California, Berkeley, approved the study and the subject gave informed consent.

%% References.

%Dellavale, D., Bonini, F., Pizzo, F., et al. Spontaneous fast-ultradian dynamics
%of polymorphic interictal events in drug-resistant focal epilepsy, Epilepsia (submitted 2023).
%Preprint freely available at DOI: https://doi.org/10.1101/2023.04.05.23288085

%Paper describing the iEEG dataset:
%Stolk, A., Griffin, S., van der Meij, R. et al. Integrated analysis of anatomical
%and electrophysiological human intracranial data. Nat Protoc 13, 1699–1723 (2018).
%DOI: 10.1038/s41596-018-0009-6
%URL: https://www.nature.com/articles/s41596-018-0009-6

%Link to the iEEG dataset:
%https://zenodo.org/record/1201560#.ZCxavXbP2Ul

%Link to the Fieldtrip tutorial using the iEEG dataset:
%https://www.fieldtriptoolbox.org/tutorial/human_ecog/

%% Changes from previous versions.

%None.

%% Tree of dependencies.

%test_node_v1.m
%
% function_zscore_v1.m
%
% function_node_v1.m
%  function_FDF_v1.m
%   function_window_v1.m
%  function_zscore_v1.m
%  function_LFDR_v1.m
%
% function_saveFig_v1.m

%% Settings for the Parallel Computing Toolbox.

%Matlab -> Parallel (in the Workspace GUI menu) -> Manage configurations ->
%select 'local' -> properties -> set the number of workers (ClusterSize) 
%equal to the number of the processor cores.

%Test if the Parallel Computing Toolbox is licensed in this computer.
%license('test','Distrib_Computing_Toolbox')

%Build a cluster using the default cluster profile and return a cluster object.
%defaultCluster = parcluster();

%Define the pool of workers (labs).
%parpool/matlabpool open
%parpool Threads
%parpool Processes

%Close or release your pool of workers (labs).
%parpool/matlabpool close
%delete(gcp('nocreate'))

%parpool/matlabpool size, or
%parpool/matlabpool('size'), or
%gcp('nocreate'),
%returns the size of the worker pool if it is open, or 0 if the pool is closed.

%Avoid "parfor" from launching a parallel pool automatically.
%parallelSettings = parallel.Settings;
%parallelSettings.Pool.AutoCreate = false;
    
%Refs:
%Matlab -> Help -> parpool
%Matlab -> Help -> matlabpool

%% Initial settings.

clear global        %Clear all global variables.
clear all           %Clear all local variables.
close all force     %Close all figures (including the wvtool figures).
format long         %Show numbers in long format.
clc                 %Clear screen.

restoredefaultpath
%pathtool

try %MATLAB 2018 and later:
    set(groot,'defaultLegendAutoUpdate','off'); %Turn off the legend auto updates.
    %For a legend in a single figure: hAxes = gca; hAxes.Legend.AutoUpdate = 'off'; 
    %Refs:
    %https://www.mathworks.com/matlabcentral/answers/320604-how-do-i-prevent-the-legend-from-auto-updating
    %https://www.mathworks.com/help/matlab/creating_plots/default-property-values.html
catch err
    disp(['The property ''defaultLegendAutoUpdate'' is not recognized in this Matlab version.']);
end %try

%Set default font:
%get(0,'DefaultAxesFontName')
%set(0,'DefaultAxesFontName','Times New Roman')

%Turn on all the warnings.
warning('on','all');

% %Turn off all the warnings.
% warning('All the warnings will be turned off.');
% warning('off','all');

%% Access path to the functions.

PATH.functions = {...
'./functions/fieldtrip/forward',...
'./functions/fieldtrip/utilities',...
'./functions/fieldtrip',...
'./functions/node',...
'./functions/lfdr',...
'./functions/windowing',...
'./functions/filtering',...
'./functions/normalization',...
'./functions/figure',...
};

addpath(PATH.functions{:}, '-begin');

%% Parameters for the data set.

%---
%Define the parameters for the datafiles.

DATASET.path = '../data';
DATASET.dataFilename = 'SubjectUCI29_data.mat';
DATASET.hdrFilename = 'SubjectUCI29_hdr.mat';
%---

%---
%Define the parameters for channels pertaining to the linear depth electrodes.

DATASET.depths = {'RAM*', 'RHH*', 'RTH*', 'ROC*', 'LAM*', 'LHH*', 'LTH*'};
DATASET.channels2show = 'RHH';
DATASET.Nchannels2show = 3;
%---

%---
%Define the index for the trial to process.

DATASET.indTrial = 18; %1 to 26.
%---

%% Parameters of band-pass filtering for detecting anomalies (i.e. amplitude outliers).

%---
%Fraction of the signal to be reflected before filtering.
NODE.reflectFraction = 1.0;
%---

%---
%Compute the cutoff frequencies of the band-pass filter for HIGH_DELTA-THETA band.
f1  = 1; %[Hz] Low cutoff frequency at -Inf dB (null amplitude).
f1p = 3; %[Hz] Low cutoff frequency at 0 dB.
f2p = 8; %[Hz] High cutoff frequency at 0 dB.
f2  = f2p + f1p - f1; %[Hz] High cutoff frequency at -Inf dB (null amplitude).  
tukey_r = 2*(f1-f1p) / (2*f1-f1p-f2p); %Parameter for the Tukey window.
%
%Frequency Domain Band-Pass Filter parameters for function "function_FDF_v1.m"
NODE.BPFcfg{1} = struct(...
'f1',f1,... %[Hz] (null amplitude = -Inf dB)
'f2',f2,... %[Hz] (null amplitude = -Inf dB)
'type','bandpass',...
'zeropadding',0,...
'freqWindowParam',struct('name','tukey','r',tukey_r),...
'timeWindowParam',struct('name','hann','sflag','symmetric'),...
'conv','circular',...
'causal',false,...
'fs',[],...
'function','function_FDF');
%
%Compute the cutoff frequencies of the band-pass filter for ALPHA-BETA band.
f1  = 8; %[Hz] Low cutoff frequency at -Inf dB (null amplitude).
f1p = 10; %[Hz] Low cutoff frequency at 0 dB.
f2p = 30; %[Hz] High cutoff frequency at 0 dB.
f2  = f2p + f1p - f1; %[Hz] High cutoff frequency at -Inf dB (null amplitude).  
tukey_r = 2*(f1-f1p) / (2*f1-f1p-f2p); %Parameter for the Tukey window.
%
%Frequency Domain Band-Pass Filter parameters for function "function_FDF_v1.m"
NODE.BPFcfg{2} = struct(...
'f1',f1,... %[Hz] (null amplitude = -Inf dB)
'f2',f2,... %[Hz] (null amplitude = -Inf dB)
'type','bandpass',...
'zeropadding',0,...
'freqWindowParam',struct('name','tukey','r',tukey_r),...
'timeWindowParam',struct('name','hann','sflag','symmetric'),...
'conv','circular',...
'causal',false,...
'fs',[],...
'function','function_FDF');
%
%Compute the cutoff frequencies of the band-pass filter for GAMMA-HIGH_GAMMA band.
f1  = 30; %[Hz] Low cutoff frequency at -Inf dB (null amplitude).
f1p = 35; %[Hz] Low cutoff frequency at 0 dB.
f2p = 150; %[Hz] High cutoff frequency at 0 dB.
f2  = f2p + f1p - f1; %[Hz] High cutoff frequency at -Inf dB (null amplitude).  
tukey_r = 2*(f1-f1p) / (2*f1-f1p-f2p); %Parameter for the Tukey window.
%
%Frequency Domain Band-Pass Filter parameters for function "function_FDF_v1.m"
NODE.BPFcfg{3} = struct(...
'f1',f1,... %[Hz] (null amplitude = -Inf dB)
'f2',f2,... %[Hz] (null amplitude = -Inf dB)
'type','bandpass',...
'zeropadding',0,...
'freqWindowParam',struct('name','tukey','r',tukey_r),...
'timeWindowParam',struct('name','hann','sflag','symmetric'),...
'conv','circular',...
'causal',false,...
'fs',[],...
'function','function_FDF');
%
%Compute the cutoff frequencies of the band-pass filter for RIPPLES.
f1  = 150; %[Hz] Low cutoff frequency at -Inf dB (null amplitude).
f1p = 155; %[Hz] Low cutoff frequency at 0 dB.
f2p = 250; %[Hz] High cutoff frequency at 0 dB.
f2  = f2p + f1p - f1; %[Hz] High cutoff frequency at -Inf dB (null amplitude).  
tukey_r = 2*(f1-f1p) / (2*f1-f1p-f2p); %Parameter for the Tukey window.
%
%Frequency Domain Band-Pass Filter parameters for function "function_FDF_v1.m"
NODE.BPFcfg{4} = struct(...
'f1',f1,... %[Hz] (null amplitude = -Inf dB)
'f2',f2,... %[Hz] (null amplitude = -Inf dB)
'type','bandpass',...
'zeropadding',0,...
'freqWindowParam',struct('name','tukey','r',tukey_r),...
'timeWindowParam',struct('name','hann','sflag','symmetric'),...
'conv','circular',...
'causal',false,...
'fs',[],...
'function','function_FDF');
%---

%% Parameters to compute the Local False Discovery Rate (LFDR) method.

%Parameters to compute the HISTOGRAM.
LFDR.histBinMethod = 'auto';
LFDR.histNormalization = 'pdf';

%Parameters to compute the MIXTURE probability density function.
LFDR.mixtureDistFunction = 'fitdist';
LFDR.mixtureDistName = 'Kernel';
LFDR.mixtureDistParam = {'Kernel', 'epanechnikov'};

%Parameters to compute the probability density function under the THEORETICAL null hypothesis.
LFDR.theoreticalDistFunction = 'makedist';
LFDR.theoreticalDistName = 'Normal';
LFDR.theoreticalDistParam = {'mu', 0, 'sigma', 1};

%Parameters to compute the probability density function under the EMPIRICAL null hypothesis.
LFDR.empiricalDistRange = []; %Use the default value.
%
LFDR.empiricalDistRangeMethod = 'auto';
LFDR.empiricalDistNintervals = 100;
%
LFDR.empiricalDistFunction = 'fitdist';
LFDR.empiricalDistName = 'Normal';
LFDR.empiricalDistParam = {};

%Define the maximum false discovery rate value (i.e. maximum rate of false positives).
LFDR.lfdr_threshold = [0.1, 0.5]; %Rate of false positives.

%Flag to plot the results.
LFDR.plotFlag = false;

%% Parameters to compute the events detection using the NODE algorithm.

NODE.LFDR = LFDR; %Structure for the local false discovery rate analysis.

%Length of the time window centered at each detected anomaly (i.e. amplitude outlier).
NODE.windowLength = 0.2; %[sec]
NODE.ampEnv = true(size(NODE.BPFcfg)); %Flags to compute the amplitude envelope in each frequency band.
NODE.fs = []; %[Hz] Sampling rate.

%Define the length of the time windows to merge the events co-occurring in time.
NODE.timeWindowEvent = NODE.windowLength; %[sec]

%Threshold for clustering.
NODE.clustThreshold = 0; %Maximum difference between the labels. Each label digit is in the range [0, 1].
   
%% Plot parameters.

%Get the screen size.
scrsz = get(0,'ScreenSize'); %scrsz -> [left, bottom, width, height].

figure_location = scrsz; %Figure in full screen mode.
%figure_location = [scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]; %[left, bottom, width, height]. 

%Define the colors.
rgbBLUE = [0,0,1];
rgbDARKBLUE  = [0,0,0.7];
rgbLIGHTBLUE  = [0.3,0.3,1];
rgbVERYLIGHTBLUE  = [0.7,0.7,1];
rgbRED = [1,0,0];
rgbDARKRED = [0.7,0,0];
rgbLIGHTRED = [1,0.3,0.3];
rgbVERYLIGHTRED = [1,0.7,0.7];
rgbGREEN = [0,1,0];
rgbDARKGREEN = [0,0.5,0];
rgbLIGHTGREEN = [0.3,0.7,0.3];
rgbVERYLIGHTGREEN = [0.5,0.7,0.5];
rgbMAGENTA = [1,0,1];
rgbCYAN = [0,1,1];
rgbYELLOW = [1,1,0];
rgbWHITE = [1,1,1];
rgbGRAY  = [0.5,0.5,0.5];
rgbDARKGRAY  = [0.3,0.3,0.3];
rgbLIGHTGRAY  = [0.7,0.7,0.7];
rgbVERYLIGHTGRAY  = [0.9,0.9,0.9];
rgbBLACK = [0,0,0];

lineProps = {...
{'Marker','none','MarkerSize',10,'MarkerEdgeColor',rgbBLUE,'MarkerFaceColor',rgbBLUE,...
'LineStyle','-','LineWidth',1,'Color',rgbBLUE},...
{'Marker','none','MarkerSize',10,'MarkerEdgeColor',rgbRED,'MarkerFaceColor',rgbRED,...
'LineStyle','-','LineWidth',1,'Color',rgbRED},...
{'Marker','none','MarkerSize',10,'MarkerEdgeColor',rgbGRAY,'MarkerFaceColor',rgbGRAY,...
'LineStyle','-','LineWidth',1,'Color',rgbGRAY},...
{'Marker','none','MarkerSize',10,'MarkerEdgeColor',rgbDARKBLUE,'MarkerFaceColor',rgbDARKBLUE,...
'LineStyle','-','LineWidth',1,'Color',rgbDARKBLUE},...
{'Marker','none','MarkerSize',10,'MarkerEdgeColor',rgbLIGHTGRAY,'MarkerFaceColor',rgbLIGHTGRAY,...
'LineStyle','-','LineWidth',1,'Color',rgbLIGHTGRAY}};

%Font size.
FONT_SIZE = 18;

%Refs:
%https://www.mathworks.com/help/matlab/ref/colororder.html

%% Define global parameters.

dimRows = 1;
dimCols = 2;

DATASET.dimChannels = 1; %rows => channels.
DATASET.dimSamples = 2; %columns => samples.

%% Set the output path and filename.

OUTPUT.fileName = ['test_node_v1'];

OUTPUT.path = ['./', OUTPUT.fileName, '/'];

if ~exist(OUTPUT.path, 'dir')
    mkdir(OUTPUT.path);
end %if ~exist(OUTPUT.path, 'dir')

%% Start the log file.

diary(strcat(OUTPUT.path, OUTPUT.fileName, '_trial', num2str(DATASET.indTrial), '.log'))  

%% Load the data.

%Load the header.
DATASET.hdr = load([DATASET.path, '/', DATASET.hdrFilename]);
DATASET.hdr = DATASET.hdr.hdr;

%Load the data.
DATASET.data = load([DATASET.path, '/', DATASET.dataFilename]);
DATASET.data = DATASET.data.data;

%% Update the configuration parameters.

%Sampling rate.
NODE.fs = DATASET.data.fsample; %[Hz]
%
for bb=1:+1:length(NODE.BPFcfg) %Loop over the frequency bands.
    NODE.BPFcfg{bb}.fs = DATASET.data.fsample; %[Hz]
end %Loop over the frequency bands.

%% Check if the sampling rate satisfy the Nyquist criterion.

for bb=1:+1:length(NODE.BPFcfg) %Loop over the frequency bands.
    assert(NODE.fs > (2*NODE.BPFcfg{bb}.f2),...
           ['The sampling rate does not meet the Nyquist criterion for frequency band ',num2str(bb),'.']);
end %Loop over the frequency bands.

%% Check the consistency between the frequency bands and the lenght of the time series.

%Compute the length of the time series.
timeSeriesLength = DATASET.data.time{DATASET.indTrial}(end) - DATASET.data.time{DATASET.indTrial}(1); %[sec]

for bb=1:+1:length(NODE.BPFcfg) %Loop over the frequency bands.
    assert(timeSeriesLength > (1/NODE.BPFcfg{bb}.f1),...
           ['The time series is shorter than one period of the frequency band ',num2str(bb),'.']);
end %Loop over the frequency bands.

%% Compute the referencing montage.

%Apply a bipolar montage to the depth electrodes.
for dd = 1:numel(DATASET.depths)
cfg            = [];
cfg.channel    = ft_channelselection(DATASET.depths{dd}, DATASET.data.label);
cfg.reref      = 'yes';
cfg.refchannel = 'all';
cfg.refmethod  = 'bipolar';
cfg.updatesens = 'yes';
DATASET.reref_depths{dd} = ft_preprocessing(cfg, DATASET.data);
end %for dd = 1:numel(depths)

%% Reset the matlab path.
%NOTE: We include this to avoid conflicts with the fieldtrip functions.

restoredefaultpath

%Add the path to the functions.
addpath(PATH.functions{:}, '-begin');

%% Extract the trial and depth electrode of interest.

for dd = 1:numel(DATASET.reref_depths)
    
    if all(contains(DATASET.reref_depths{dd}.label,DATASET.channels2show,'IgnoreCase',true))

        %Extract the time vector (1 x samples).
        DATASET.t = DATASET.reref_depths{dd}.time{DATASET.indTrial};

        %Extract the signal (channels x samples).
        DATASET.signal = DATASET.reref_depths{dd}.trial{DATASET.indTrial};
        
        %Extract the labels of the channels.
        DATASET.channelLabel = DATASET.reref_depths{dd}.label;

        break; %Stop the loop.

    end %if...

end %for dd = 1:numel(DATASET.reref_depths)

%Compute the number of samples.
DATASET.Nsamples = size(DATASET.signal,DATASET.dimSamples); %[samples]

%% Z-score normalization of the current epoch.

%Z-score normalization of each channel (rows) across the samples (columns).
DATASET.signal = function_zscore_v1(DATASET.signal.').';

%% Loop over the groups of channels.

%IMPORTANT:
%In this script we process all the channels taken individually. That is, only
%one channel is included in each group of channels. 

%Initialize the markers and the kernels.
NODE.results.markers = [];
%
for chgroup=1:+1:size(DATASET.signal,DATASET.dimChannels) %Loop over the groups of channels.
    
%% Update the output filename.

OUTPUT.fileName_chgroup = [DATASET.channelLabel{chgroup}, '_trial', num2str(DATASET.indTrial)];
    
%% Extract the channels of interest.

%Display the channel being processed.
disp(['Processing channel ', DATASET.channelLabel{chgroup}, '; ', num2str(chgroup),' out of ',num2str(size(DATASET.signal,DATASET.dimChannels)),'.']);

%Extract a single channel.
DATASET.signal_chgroup = DATASET.signal(chgroup,:);

%Label of the channels.
NODE.channelLabel = DATASET.channelLabel(chgroup);

%% Detection of the events using the NODE algorithm.

%[Dout] = function_node_v1(DATASET.signal_chgroup, NODE);
%[Dout, Sout] = function_node_v1(DATASET.signal_chgroup, NODE);
[Dout, Sout, BPFout] = function_node_v1(DATASET.signal_chgroup, NODE);

%---
%How to read the output structures:
%cc -> index for channels.
%bb -> index for the frequency bands.
%th -> index for LFDR thresholds.
%ee -> index for events.
%aa -> index for anomalies.
%
%Dout.ampThreshold(cc,bb,th)
%Dout.markers.anomalyPosition{cc,bb,th}(aa)
%Dout.markers.eventPosition{cc,th}(ee)
%Dout.markers.eventLabel{cc,th}{ee}
%Dout.Dcfg
%
%Sout.eventSignal{cc,th}{ee}(samples)
%Sout.eventT{cc,th}{ee}(samples)
%
%BPFout.eventBPFsignal{cc,th}{ee}(samples,bb)
%BPFout.eventT{cc,th}{ee}(samples)
%BPFout.eventF{bb}
%---

%% Loop over the LFDR thresholds.

for th=1:+1:length(LFDR.lfdr_threshold) %Loop over the LFDR thresholds.

%Display the number of the local false discovery rate being processed.
disp(['Processing local FDR ',num2str(th),' out of ',num2str(length(LFDR.lfdr_threshold)),'.']);    

%% Add the time offset to the events corresponding to the current epoch.

for cc=1:+1:size(Dout.markers.eventPosition,dimRows) %Loop over the channels.
    
    Dout.markers.eventPosition{cc,th} = Dout.markers.eventPosition{cc,th} + DATASET.t(1);
    %
    for bb=1:+1:length(NODE.BPFcfg) %Loop over the frequency bands.
        Dout.markers.anomalyPosition{cc,bb,th} = Dout.markers.anomalyPosition{cc,bb,th} + DATASET.t(1);
    end %Loop over the frequency bands.    

    for ee=1:+1:length(Dout.markers.eventLabel{cc,th}) %Loop over the events.    
        
        if exist('Sout','var') && ~isempty(Sout.eventT{cc,th}{ee})
            Sout.eventT{cc,th}{ee} = Sout.eventT{cc,th}{ee} + DATASET.t(1);
        end %if exist('Sout','var') && ~isempty(Sout.eventT{cc,th}{ee})
        
        if exist('BPFout','var') && ~isempty(BPFout.eventT{cc,th}{ee})
            BPFout.eventT{cc,th}{ee} = BPFout.eventT{cc,th}{ee} + DATASET.t(1);
        end %if exist('BPFout','var') && ~isempty(BPFout.eventT{cc,th}{ee})
        
    end %Loop over the events.

end %Loop over the channels. 

%% Save the amplitude threshold corresponding to the current LFDR threshold.

%Dout.ampThreshold         -> (matrix: channels x frequency bands x LFDR thresholds)
%NODE.results.ampThreshold -> (matrix: channels x frequency bands x LFDR thresholds)

NODE.results.ampThreshold(:,:,th) = Dout.ampThreshold(:,:,th);

%% Save the markers of the current LFDR threshold.

if (th == 1)
    
%Initialize the markers.
NODE.results.markers(cc).eventPosition = Dout.markers.eventPosition{cc,th};
NODE.results.markers(cc).eventLabel = Dout.markers.eventLabel{cc,th};
%
for cc=1:+1:size(Dout.markers.eventPosition,dimRows) %Loop over the channels.
    NODE.results.markers(cc).indThreshold = th * ones(size(Dout.markers.eventLabel{cc,th}));
end %Loop over the channels.

else %if (th == 1)
    
%Update the markers with the current events.
%
for cc=1:+1:size(Dout.markers.eventPosition,dimRows) %Loop over the channels.
    NODE.results.markers(cc).eventPosition = [NODE.results.markers(cc).eventPosition, Dout.markers.eventPosition{cc,th}];
    NODE.results.markers(cc).eventLabel = [NODE.results.markers(cc).eventLabel, Dout.markers.eventLabel{cc,th}];
    NODE.results.markers(cc).indThreshold = [NODE.results.markers(cc).indThreshold, th * ones(size(Dout.markers.eventLabel{cc,th}))];
end %Loop over the channels.

end %if (th == 1)

end %Loop over the LFDR thresholds.

%% Merge the events co-occurring in time.

%---
%IMPORTANT: Note on the clustering based on iterations over the local LFDR thresholds.
% 
%It is important to note that we could implement the consistency clustering
%using the full resolution fingerprint "logicalLabel" corresponding to a given event: 
% 
%logicalLabel (Num. of local FDR values x Num. of frequency bands = 2 x 4 logical array):
%0100
%1110
%
%A clustering using the matching between the "logicalLabel" of each event will 
%result in a more fine classification of the events with respect to use the 
%weighted average fingerprint "eventDecimalLabel".
%---

%Sort the events according to ascending position.
for cc=1:+1:length(NODE.results.markers) %Loop over the channels.
    [NODE.results.markers(cc).eventPosition, indSort] = sort(NODE.results.markers(cc).eventPosition,'ascend');
    NODE.results.markers(cc).eventLabel = NODE.results.markers(cc).eventLabel(indSort);
    NODE.results.markers(cc).indThreshold = NODE.results.markers(cc).indThreshold(indSort);
end %Loop over the channels.
%
%Refs:
%https://fr.mathworks.com/help/matlab/ref/double.cat.html
%https://fr.mathworks.com/help/matlab/ref/sort.html

%Compute the label and consistency for each event.
for cc=1:+1:length(NODE.results.markers) %Loop over the channels.
    
    %Initialize the variables.
    position = NODE.results.markers(cc).eventPosition;
    label = NODE.results.markers(cc).eventLabel;
    indThreshold = NODE.results.markers(cc).indThreshold;
    %
    NODE.results.markers(cc).eventPosition = [];
    NODE.results.markers(cc).decimalLabelEvents = [];
    if isempty(position)
        %To assign this default value is necessary since it is not possible to
        %concatenate empty 3D-matrices (Nrows x Ncolumns x Nslices).
        NODE.results.markers(cc).binaryLabelEvents = NaN(length(LFDR.lfdr_threshold),length(NODE.BPFcfg));
    else
        NODE.results.markers(cc).binaryLabelEvents = [];
    end %if isempty(position)
    %    
    ee = 1; %Initialize the index for the events.
    while ~isempty(position) %Loop over the events.
        
        %Check if any DETECTION is within the time window of the current CURRENT DETECTION.
        indStartEvent = position > position(ee)-NODE.timeWindowEvent/2;
        indEndEvent = position < position(ee)+NODE.timeWindowEvent/2;
        %
        indNewEvent = indStartEvent & indEndEvent;
        
        %Merge the co-occurring events in a new single event.
        NODE.results.markers(cc).eventPosition = [NODE.results.markers(cc).eventPosition; mean(position(indNewEvent))];
        %
        %Initialize the variables.
        charLabel = char(label{indNewEvent});
        doubleLabel = NaN(1,size(charLabel,dimCols));
        logicalLabel = false(length(LFDR.lfdr_threshold),length(NODE.BPFcfg));
        %Compute the labels.        
        for digit=1:+1:size(charLabel,dimCols)
            %
            %Compute the logical label.
            logicalLabel(indThreshold(indNewEvent),digit) = logical(str2num(charLabel(:,digit)));
            %logicalLabel(indThreshold(indNewEvent),digit) = logical(charLabel(:,digit) - '0');
            %Refs: 
            %https://fr.mathworks.com/matlabcentral/answers/112717-convert-class-char-to-class-logical-how
            %https://stackoverflow.com/questions/20032413/matlab-string-vector-character-subtraction
            %
            %Compute the upper-bound label.
            doubleLabel(digit) = max(logicalLabel(:,digit) .* (1-LFDR.lfdr_threshold.'));
            %Compute the weighted average label.
            %doubleLabel(digit) = sum(logicalLabel(:,digit) .* (1-LFDR.lfdr_threshold.'));
            %Compute the average label.
            %doubleLabel(digit) = sum(logicalLabel(:,digit));
            %
        end %for digit=1:+1:size(charLabel,dimCols)
        %Compute the decimal label based on the upper-bound label.
        NODE.results.markers(cc).decimalLabelEvents = [NODE.results.markers(cc).decimalLabelEvents; doubleLabel];        
        %Compute the decimal label based on the weighted average label.
        %NODE.results.markers(cc).decimalLabelEvents = [NODE.results.markers(cc).decimalLabelEvents; doubleLabel/sum(1-LFDR.lfdr_threshold)];
        %Compute the decimal label based on the average label.
        %NODE.results.markers(cc).decimalLabelEvents = [NODE.results.markers(cc).decimalLabelEvents; doubleLabel/length(LFDR.lfdr_threshold)];
        %
        %Compute the binary label.
        NODE.results.markers(cc).binaryLabelEvents = cat(3, NODE.results.markers(cc).binaryLabelEvents, logicalLabel);
        %Refs:
        %https://fr.mathworks.com/help/matlab/ref/double.cat.html
        
        %Clear the events that were merged in a new single event,
        %so the next event is always the first event -> ee = 1;
        position(indNewEvent) = [];
        label(indNewEvent) = [];
        indThreshold(indNewEvent) = [];
        
        %The next event is always the first event.
        ee = 1;        

    end %Loop over the events.
end %Loop over the channels.
%
%Remove the originals labels from the markers structure.
%
if isfield(NODE.results.markers,'label')
    NODE.results.markers = rmfield(NODE.results.markers,'label');
end %if isfield(NODE.results.markers,'label')
%
if isfield(NODE.results.markers,'indThreshold')
    NODE.results.markers = rmfield(NODE.results.markers,'indThreshold');
end %if isfield(NODE.results.markers,'indThreshold')

%% Unwrap the markers over the channels.

%NOTE: Since here the events pertaining to different channels are
%merged in a single array, the difference between succesive positions is
%not longer determined by the parameter NODE.windowLength.

%Initialize the variables.
NODE.results.eventPosition = [];
NODE.results.channel = [];
NODE.results.decimalLabelEvents = [];
NODE.results.binaryLabelEvents = [];
for cc=1:+1:length(NODE.results.markers) %Loop over the channels.
    
    %Register the channel number for the current events.
    NODE.results.channel = [NODE.results.channel; cc*ones(length(NODE.results.markers(cc).eventPosition),1)];
    %Unwrap the position of the detected events [sec].
    NODE.results.eventPosition = [NODE.results.eventPosition; NODE.results.markers(cc).eventPosition];
    %Unwrap the decimal label of the detected events.
    NODE.results.decimalLabelEvents = [NODE.results.decimalLabelEvents; NODE.results.markers(cc).decimalLabelEvents];
    %Unwrap the binary label of the detected events.
    NODE.results.binaryLabelEvents = cat(3, NODE.results.binaryLabelEvents, NODE.results.markers(cc).binaryLabelEvents);

end %Loop over the channels.

%Sort the events according to ascending position.
[NODE.results.eventPosition, indSort] = sort(NODE.results.eventPosition,'ascend');
NODE.results.channel = NODE.results.channel(indSort);
NODE.results.decimalLabelEvents = NODE.results.decimalLabelEvents(indSort,:);
NODE.results.binaryLabelEvents = NODE.results.binaryLabelEvents(:,:,indSort);
%
%Refs:
%https://fr.mathworks.com/help/matlab/ref/double.cat.html
%https://fr.mathworks.com/help/matlab/ref/sort.html

%% Clear the field containing the wrapped markers.

if isfield(NODE.results,'markers')
    NODE.results = rmfield(NODE.results,'markers');
end %if isfield(NODE.results,'markers')

%% Remove the events likely associated with edge artifacts.

indEvents = ( NODE.results.eventPosition > (DATASET.t(1) + NODE.timeWindowEvent/2) ) &...
            ( NODE.results.eventPosition < (DATASET.t(end) - NODE.timeWindowEvent/2) );

NODE.results.eventPosition = NODE.results.eventPosition(indEvents);
NODE.results.channel = NODE.results.channel(indEvents);
NODE.results.decimalLabelEvents = NODE.results.decimalLabelEvents(indEvents,:);
NODE.results.binaryLabelEvents = NODE.results.binaryLabelEvents(:,:,indEvents);

%% Compute the cluster indices using the event labels.

%Initialize the variables.
NODE.results.idEvents = [];
NODE.results.decimalLabelClusters = [];
NODE.results.Nec = [];
%  
decimalLabelEvents = NODE.results.decimalLabelEvents;
%
NODE.results.idEvents = NaN(size(decimalLabelEvents,dimRows),1);
uu = 1; ee = 1; %Initialize the indices for the clusters and the events.
while ~isempty(decimalLabelEvents) %Loop over the events.
    
    %Compute the indices for the events with the same label.
    %indDetections = fix(sum(NODE.results.decimalLabelEvents == decimalLabelEvents(ee,:),dimCols)/size(decimalLabelEvents,dimCols)) == 1;
    indDetections = all( abs(NODE.results.decimalLabelEvents - decimalLabelEvents(ee,:)) <= NODE.clustThreshold, dimCols);    
    
    %Merge those events in the same cluster.
    NODE.results.idEvents(indDetections) = uu;
    
    %Compute the decimal label for the current cluster.
    NODE.results.decimalLabelClusters(uu,:) = decimalLabelEvents(ee,:);
    %Compute the number of events in the current cluster. 
    NODE.results.Nec(uu,1) = sum(indDetections);
    
    %Increase the cluster index.
    uu = uu + 1;
        
    %Clear the events that were merged in a single cluster,
    %so the next event is always the first event -> ee = 1;
    %indKeep = not( fix(sum(decimalLabelEvents == decimalLabelEvents(ee,:),dimCols)/size(decimalLabelEvents,dimCols)) == 1 );
    indKeep = not( all( abs(decimalLabelEvents - decimalLabelEvents(ee,:)) <= NODE.clustThreshold, dimCols) );    
    decimalLabelEvents = decimalLabelEvents(indKeep,:);
    
    %The next event is always the first event.
    ee = 1;        
    
end %Loop over the events.

%Compute the resulting number of clusters.
NODE.results.Ncluster = numel(unique(NODE.results.idEvents));   

%Compute the resulting number of detected events.
assert(sum(NODE.results.Nec)==length(NODE.results.idEvents),['Inconsistency detected while computing the number of detected events.']);
NODE.results.Nevent = sum(NODE.results.Nec);

%---
%Sort the clusters according to the descending order of number of elements.
[NODE.results.Nec, indClust] = sort(NODE.results.Nec, 'descend');
%Sort the decimal label of the clusters accordingly.
NODE.results.decimalLabelClusters = NODE.results.decimalLabelClusters(indClust,:);
%Recompute the cluster indices accordingly.
for uu=1:+1:NODE.results.Ncluster
    %Remplace the value of the cluster indices using an offset,
    %so the new sorted indices are different from the original ones. 
    NODE.results.idEvents(NODE.results.idEvents==indClust(uu)) = uu + NODE.results.Ncluster;
end %for uu=1:+1:NODE.results.Ncluster
%Remove the offset.
NODE.results.idEvents = NODE.results.idEvents - NODE.results.Ncluster;
%---

%---
%Sort the clusters according to the non-null label digit in the same time scales.

aux_decimalLabelClusters = [];
aux_Nec = [];
aux_indClust = [];
aux_idClusters = []; indCoC = 1;
for uu=1:+1:NODE.results.Ncluster %Loop over the clusters.

if any(NODE.results.decimalLabelClusters(uu,:)) %Only process not all-null label clusters.

%Compute the indices for the clusters with non-null label digit in the same time scales.    
indSubClust = all(not(xor(NODE.results.decimalLabelClusters(uu,:),NODE.results.decimalLabelClusters)),dimCols);

%Compute the decimal label. 
aux_decimalLabelClusters = [aux_decimalLabelClusters; NODE.results.decimalLabelClusters(indSubClust,:)];
%Compute the number of events in each cluster.
aux_Nec = [aux_Nec; NODE.results.Nec(indSubClust)];
%Compute the indices for the current clusters.
aux_indClust = [aux_indClust; find(indSubClust)];

%Compute the label for the cluster of clusters (CoC), 
%i.e: Cluster with non-null label digit in the same time scales.
aux_idClusters = [aux_idClusters; indCoC*ones(sum(indSubClust),1)];
indCoC = indCoC + 1;

%Assign a null in all time scales to avoid re-processing the current clusters.
NODE.results.decimalLabelClusters(indSubClust,:) = 0;

end %if any(NODE.results.decimalLabelClusters(uu,:))

end %Loop over the clusters.

%Assign the re-ordered clusters.
NODE.results.decimalLabelClusters = aux_decimalLabelClusters;
NODE.results.Nec = aux_Nec;
NODE.results.idClusters = aux_idClusters;

%Recompute the cluster indices accordingly.
for uu=1:+1:NODE.results.Ncluster %Loop over the clusters.
    %Remplace the value of the cluster indices using an offset,
    %so the new sorted indices are different from the original ones. 
    NODE.results.idEvents(NODE.results.idEvents==aux_indClust(uu)) = uu + NODE.results.Ncluster;
end %Loop over the clusters.
%Remove the offset.
NODE.results.idEvents = NODE.results.idEvents - NODE.results.Ncluster;
%---

%% Save the results.

%Save the output data using the default mat file version (v6: MATLAB 5.0 MAT-file).
%
%For variables larger than 2GB use MAT-file version 7.3 or later.
save([OUTPUT.path, OUTPUT.fileName_chgroup],...
     'NODE','-v6'); %Structure containing the information of the events.
     %'NODE','-v7.3'); %Structure containing the information of the events.
%
% if exist([OUTPUT.path, OUTPUT.fileName_chgroup, '.mat'],'file')
%     warning('MATLAB:test_node',['The output .mat file already exists. File not saved!']);
% else
%     %For variables larger than 2GB use MAT-file version 7.3 or later.
%     save([OUTPUT.path, OUTPUT.fileName_chgroup],...
%          'NODE','-v6'); %Structure containing the information of the events.
%          %'NODE','-v7.3'); %Structure containing the information of the events.
% end %if exist([OUTPUT.path, OUTPUT.fileName_chgroup, '.mat'],'file')

%Check the contents of the .mat file using the whos function.
%disp(['Contents of:',sprintf('\n'),ofile]);
%whos('-file',ofile)

%Get the info of the ".mat" file.
%type strcat(ofile, '.mat')
%
%Refs:
%Matlab -> Help -> save
%Matlab -> Help -> mat file
%http://www.mathworks.com/help/matlab/import_export/mat-file-versions.html

end %Loop over the groups of channels.

%% ========================================================================

%% Clear the variables associated with individual group of channels.

if isfield(DATASET,'signal_chgroup')
    DATASET = rmfield(DATASET,'signal_chgroup');
end %if isfield(DATASET,'signal_chgroup')

if isfield(OUTPUT,'fileName_chgroup')
    OUTPUT = rmfield(OUTPUT,'fileName_chgroup');
end %if isfield(OUTPUT,'fileName_chgroup')

NODE = [];

%% Loop over the groups of channels.

%IMPORTANT:
%In this script we process all the channels taken individually. That is, only
%one channel is included in each group of channels. 

%Define the number of channels to show.
%DATASET.Nchannels2show = size(DATASET.signal,DATASET.dimChannels);

%Initialize the variables.
NODE.results.eventPosition = {};
NODE.results.decimalLabel = {};
NODE.results.strDecimalLabel = {};
NODE.results.channel = {};
%
for chgroup=1:+1:DATASET.Nchannels2show %Loop over the groups of channels.

%% Load the .mat file with the markers.

%Define the file name.
fileName = [OUTPUT.path, DATASET.channelLabel{chgroup}, '_trial', num2str(DATASET.indTrial), '.mat'];

%Clear the data from the previous epoch.
NODE.node = [];
%
try 
    %Try to load the file for the current epoch.
    NODE.node = load([fileName], '-mat');
    NODE.node = NODE.node.NODE;
catch err
    %Display a message.
    disp([err.message]);
end %try

%% Extract some relevant parameters.   

%Sampling rate.
NODE.fs = NODE.node.fs; %[Hz]

%Length of the time window centered at each detected event.
NODE.windowLength = NODE.node.windowLength; %[sec]

%Compute the label of the channel.
NODE.results.channelLabel(chgroup) = NODE.node.channelLabel;

%% Extract the info of the detected events.

NODE.results.eventPosition{chgroup} = NODE.node.results.eventPosition;
NODE.results.decimalLabel{chgroup} = NODE.node.results.decimalLabelEvents;
NODE.results.channel{chgroup} = NODE.node.results.channel;    

%Compute the string for the decimal labels.
for ee=1:+1:length(NODE.results.eventPosition{chgroup}) %Loop over the events.

    %Synthesize the string for the decimal label.
    strDecimalLabel = [];
    for dd=1:+1:length(NODE.results.decimalLabel{chgroup}(ee,:))
        strDecimalLabel = [strDecimalLabel, sprintf('%.4g',NODE.results.decimalLabel{chgroup}(ee,dd)), '_']; 
    end %for dd=1:+1:length(NODE.results.decimalLabel{chgroup}(ee,:))
    %Remove the extra '_' character at the end of the string.
    strDecimalLabel = strDecimalLabel(1:end-1);
    %Remove the dots.
    strDecimalLabel = strrep(strDecimalLabel,'.','');
    
    %Compute the string for the decimal labels.
    NODE.results.strDecimalLabel{chgroup}{ee,1} = strDecimalLabel;
    
end %Loop over the events.

end %Loop over the groups of channels.

%% ========================================================================

%% Show the events for the clusters and channels of interest.

%Offset of the amplitude to show the channels stacked in the same plot.
%OFFSET_AMP = max(abs(DATASET.signal(:)));
OFFSET_AMP = 1.25;

%Create a new figure.
hFigure = figure;
%hFigure = figure('Position',figure_location);
hold on,

%---
%Plot the quantities of interest.
for chgroup=1:+1:DATASET.Nchannels2show %Loop over the groups of channels.

    %---
    %Plot the horizontal reference line.
    plot([DATASET.t(1), DATASET.t(end)], [-OFFSET_AMP*(chgroup-1), -OFFSET_AMP*(chgroup-1)],...
         'LineStyle','--','LineWidth',1,'Color',rgbBLACK);
    %---    
    
    %---
    %Plot the signal.
    plot(DATASET.t, DATASET.signal(chgroup,:)/max(abs(DATASET.signal(chgroup,:)))-OFFSET_AMP*(chgroup-1), lineProps{5}{:}, 'LineWidth',2);
    %---
    
    %---
    if ~isempty(NODE.results.eventPosition{chgroup})
    
    %Plot the DETECTED EVENTS.
    y = -OFFSET_AMP*(chgroup-1); ee = 1; 
    ph(1) =...
    plot(NODE.results.eventPosition{chgroup}(ee),y,'v','MarkerSize',8,'MarkerEdgeColor',rgbRED,'MarkerFaceColor',rgbRED);
    %
    %Plot the TIME WINDOW around the DETECTED EVENTS.
    plot([NODE.results.eventPosition{chgroup}(ee)-NODE.windowLength/2,...
         NODE.results.eventPosition{chgroup}(ee)+NODE.windowLength/2],...
         [y, y],'LineStyle','-','LineWidth',2,'Color',rgbRED);
    %
    plot(NODE.results.eventPosition{chgroup}(ee)-NODE.windowLength/2,y,'>','MarkerSize',8,'MarkerEdgeColor',rgbRED,'MarkerFaceColor',rgbRED);
    plot(NODE.results.eventPosition{chgroup}(ee)+NODE.windowLength/2,y,'<','MarkerSize',8,'MarkerEdgeColor',rgbRED,'MarkerFaceColor',rgbRED);
    %
    %Show the labels.
    text(NODE.results.eventPosition{chgroup}(ee),y,...
    [NODE.results.strDecimalLabel{chgroup}{ee}],...
    'HorizontalAlignment','left',...
    'VerticalAlignment','middle',...
    'Rotation',90,...
    'Color',rgbRED,...
    'fontname','Times New Roman','fontsize',FONT_SIZE,'fontweight','normal','interpreter','none');

    for ee=2:+1:length(NODE.results.eventPosition{chgroup}) %Loop over the DETECTED EVENTS.
        if NODE.results.eventPosition{chgroup}(ee) <= (NODE.results.eventPosition{chgroup}(ee-1) + NODE.windowLength/2)
            y = y + 0.1;
        else
            y = -OFFSET_AMP*(chgroup-1);
        end
        %Plot the DETECTED EVENTS.
        plot(NODE.results.eventPosition{chgroup}(ee),y,'v','MarkerSize',8,'MarkerEdgeColor',rgbRED,'MarkerFaceColor',rgbRED);
        %
        %Plot the TIME WINDOW around the DETECTED EVENTS.
        plot([NODE.results.eventPosition{chgroup}(ee)-NODE.windowLength/2,...
             NODE.results.eventPosition{chgroup}(ee)+NODE.windowLength/2],...
             [y, y],'LineStyle','-','LineWidth',2,'Color',rgbRED);
        %
        plot(NODE.results.eventPosition{chgroup}(ee)-NODE.windowLength/2,y,'>','MarkerSize',8,'MarkerEdgeColor',rgbRED,'MarkerFaceColor',rgbRED);
        plot(NODE.results.eventPosition{chgroup}(ee)+NODE.windowLength/2,y,'<','MarkerSize',8,'MarkerEdgeColor',rgbRED,'MarkerFaceColor',rgbRED);         
        %
        %Show the labels.
        text(NODE.results.eventPosition{chgroup}(ee), y,...
        [NODE.results.strDecimalLabel{chgroup}{ee}],...
        'HorizontalAlignment','left',...
        'VerticalAlignment','middle',...
        'Rotation',90,...
        'Color',rgbRED,...
        'fontname','Times New Roman','fontsize',FONT_SIZE,'fontweight','normal','interpreter','none');        
    end %Loop over the DETECTED EVENTS.
    
    end %if ~isempty(position)
    %---

end %Loop over the groups of channels.
%---

% if exist('ph','var')
% legend(ph,...
%        {'Detections',...
%         },...
%        'location','best','interpreter','none');
% %legend('boxoff')
% legend('boxon') 
% %Clear the plot handle.
% clear ph
% end %if exist(ph,'var')

box on, grid off,
axis tight;
%axis([]);
haxes = gca; 
haxes.YLim = [-OFFSET_AMP*DATASET.Nchannels2show, OFFSET_AMP];
%haxes.XLim = [];
%haxes.YLim = [];

set(gca,'xscale','linear');
set(gca,'yscale','linear');
%set(gca,'xscale','log');
%set(gca,'yscale','log');

haxes.TickLabelInterpreter = 'none';
% set(gca,'XTick',[]);
% set(gca,'XTickLabel',{'';''});
set(gca,'YTick',OFFSET_AMP*(-DATASET.Nchannels2show+1:+1:0));
set(gca,'YTickLabel',flipud(NODE.results.channelLabel(:)));

set(gca,'TickDirMode','manual','TickDir','in');  
set(gca,'YMinorTick','Off','XMinorTick','On','Linewidth',2);
set(gca,'fontname','Times New Roman','fontsize',FONT_SIZE,'fontweight','normal');

xlabel('Time [sec]','interpreter','none');
%ylabel('Channel','interpreter','none');
title([DATASET.dataFilename],'fontweight','normal','interpreter','none');
 
%Save the figure.
function_saveFig_v1(gcf, [OUTPUT.path, DATASET.channels2show, '_trial', num2str(DATASET.indTrial)], 'fig');
%function_saveFig_v1(gcf, [OUTPUT.path, DATASET.channels2show, '_trial', num2str(DATASET.indTrial)], 'png');

%% Save a copy of the running script (.cfg file).

%cfgfile = strcat(mfilename('fullpath'), '_trial', num2str(DATASET.indTrial), '.cfg');
cfgfile = strcat(OUTPUT.path, mfilename, '_trial', num2str(DATASET.indTrial), '.cfg');
%cfgfile = strcat(OUTPUT.path, OUTPUT.fileName, '.cfg');

copyfile(strcat(mfilename('fullpath'),'.m'), cfgfile);
%
% if exist(cfgfile,'file')
%     warning('MATLAB:test_node',...
%             ['The .cfg file already exists for the current script. File not saved!']);
% else
%     copyfile(strcat(mfilename('fullpath'),'.m'), cfgfile);
% end %if exist(cfgfile,'file')

%Refs:
%https://www.mathworks.com/help/matlab/ref/mfilename.html

%% End the log file.

diary off

%Refs:
%https://www.mathworks.com/help/matlab/ref/diary.html

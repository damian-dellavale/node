%%

% Departamento de Fisica Medica (DFM)
% Centro Atomico Bariloche, Comision Nacional de Energia Atomica (CNEA)
% Instituto Balseiro, Universidad Nacional de Cuyo
% Consejo Nacional de Investigaciones Científicas y Técnicas (CONICET)
% R8402AGP San Carlos de Bariloche, Rio Negro, Argentina
% https://fisica.cab.cnea.gov.ar/fisicamedica/sites-ddellavale/

% Dynamical Brain Mapping Group (DYNAMAP)
% Institut de Neurosciences des Systemes (INS)
% Faculte de Medecine, Aix-Marseille Universite
% 27, Boulevard Jean Moulin, 13005 Marseille, France
% https://ins-amu.fr/dynamap

% Authors: Damian Dellavale (dellavale@cab.cnea.gov.ar,
%                            dellavaledamian@gmail.com,
%                            hector.dellavale-clara@univ-amu.fr).

% Project: Nested Outlier Detection (NODE) algorithm.

% Date: 06/04/2023

%% Description.

%The NODE algorithm define the events by detecting anomalies (i.e. outliers)
%of amplitude across the frequency bands of interest. For this, the NODE 
%algorithm uses the Local False Discovery Rate (LFDR) method to identify 
%the transient excursions of amplitudes corresponding to outliers, while 
%controlling for the proportion of false positives. Then, the anomalies 
%co-ocurring in time are merged across the frequency bands to conform a 
%single event. The time resolution used to define co-occurence in time is 
%given by the time window args.windowLength (array: 1 x frequency bands).

%% References.

%Dellavale, D., Bonini, F., Pizzo, F., et al. Spontaneous fast-ultradian dynamics
%of polymorphic interictal events in drug-resistant focal epilepsy, Epilepsia (submitted 2023).
%Preprint freely available at DOI: https://doi.org/10.1101/2023.04.05.23288085

%% Changes from previous versions.

%None.

%% How to test this script.

%Find the instructions in the link:
%
%https://github.com/damian-dellavale/node/

%% Tree of dependencies.

%function_node_v1.m
%
% function_FDF_v1.m
%  function_window_v1.m
%
% function_zscore_v1.m
%
% function_LFDR_v1.m

%% Main function.

function [Dout, Sout, BPFout] = function_node_v1(signal, varargin)
%==========================================================================
%Inputs:
%
%signal -> Input signal (matrix: channels x samples = Nchannels x Nsamples).
%
%Dcfg   -> Configuration structure with fields:
%
%
%NOTE: 
%The input arguments can be passed in a single structure, i.e. function_node_v1(signal, Dcfg),
%or in the name-value pair format, i.e. function_node_v1(signal, 'fs', 1000,...). 

%Outputs:
%
%Dout   -> Output structure with fields:
%
%Sout   -> Output structure with fields:
%
%BPFout -> Output structure with fields:
%
%NOTE: 
%The output structures Sout and BPFout containing the time series corresponding
%to each detected event, respectively, are computed only if they are asked 
%as outputs in the current call to the function.
%==========================================================================

%Check the input arguments ------------------------------------------------

args = checkInputs(signal, varargin);

%Clear the duplicated variables outside the args structure.
signal = [];

%--------------------------------------------------------------------------

%Initialize parameters and settings ---------------------------------------

dimRows = 1;
dimCols = 2;

dimChannels = 1; %rows => channels.
dimSamples = 2; %columns => samples.

%Compute the number of channels.
args.Nchannels = size(args.signal,dimChannels);

%Compute the number of samples of the input time series.
args.Nsamples = size(args.signal,dimSamples);

%Compute the length of the input time series.
args.signalLength_sec = args.Nsamples/args.fs; %[sec]

%Compute the number of frequency bands.
NfreqBands = length(args.BPFcfg);

%Compute the number of FDR threshold values.
Nthresholds = length(args.LFDR.lfdr_threshold);

%Use the saveSignal and saveTF flag to record the time series and 
%time-frequency maps corresponding to each detected event, only
%if they were asked as outputs in the current call to the function.
%
%NOTE: The use of nargout is not allowed inside the parfor loop.
saveSignal = false;
saveBPFsignal = false;
if nargout == 3
    saveSignal = true;
    saveBPFsignal = true;
elseif nargout == 2
    saveSignal = true;
end %if nargout == 3
%Refs:
%https://fr.mathworks.com/help/matlab/ref/nargout.html
%
%https://fr.mathworks.com/matlabcentral/answers/353759-check-which-output-arguments-are-requested-within-a-function
%https://fr.mathworks.com/matlabcentral/answers/476440-how-to-check-if-one-of-output-variables-is-not-called
%https://fr.mathworks.com/matlabcentral/answers/379731-how-to-check-if-output-argument-is-ignored
%https://fr.mathworks.com/matlabcentral/answers/184617-how-to-detect-in-output-argument-list-so-as-to-avoid-memory-allocation-inside-called-function
%https://fr.mathworks.com/matlabcentral/answers/142733-equivalent-of-inputname-for-output-variables

%--------------------------------------------------------------------------

%Default values of the outputs --------------------------------------------
%--------------------------------------------------------------------------

%Configure the parallel computing toolbox ---------------------------------

%If the parallel computing toolbox is licensed and installed in this computer.
if license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))

    if ~isempty(gcp('nocreate')) %If a pool exist.

        localPool = gcp();
        Ncores = localPool.NumWorkers;

    else %If no pool exist.

        %Build a cluster using the default cluster profile and return a cluster object. 
        defaultCluster = parcluster();
           
        %Define the number of cores.
        %
%         if defaultCluster.NumWorkers >= 4
%             Ncores = min(defaultCluster.NumWorkers-1,args.Nchannels);
%         else
%             Ncores = 2;
%         end %if defaultCluster.NumWorkers >= 4
        %
        Ncores = defaultCluster.NumWorkers;

        %Create a parallel pool of workers on a cluster and return a pool object.
        localPool = parpool('local',Ncores);

    end %if ~isempty(gcp('nocreate')) %If a pool exist. 

else %If the parallel computing toolbox is not licensed in this computer.

    Ncores = 1; %Use a single core.

end %if license('test','Distrib_Computing_Toolbox')...

%--------------------------------------------------------------------------

%Configure the RAM memory utilization -------------------------------------

if ismac %Code to run on Mac platform.
    error('MATLAB:Platform', ['Please, add the code to obtain the available memory in MAC platform.']);
elseif isunix %Code to run on Linux platform.
    [~, txt] = system('top -l 1 -s 0 | grep PhysMem: | awk ''{print $6}''');
    res = str2double(strrep(txt, 'M', ''));
    RAMpCore_MB = res / Ncores;
elseif ispc %Code to run on Windows platform.
    [~, RAM] = memory;
    RAMpCore_MB = (RAM.PhysicalMemory.Available/1e6) / Ncores;
else %Unknown.
    error('MATLAB:Platform', ['OS platform not recognized.']);
end
%Ref: https://www.mathworks.com/help/matlab/ref/isunix.html

%Compute the memory required for Double-Precision Floating Point (8 bytes)
%and complex (x 2) data.
RAMrequired_MB = args.signalLength_sec * args.fs * 8 * 2 / 1e6;

%--------------------------------------------------------------------------

%Start analysis -----------------------------------------------------------

%---
%Compute the number of chunks.
Nchunks = ceil(RAMrequired_MB / RAMpCore_MB);
%---

%---
%Configure the start and duration (in sec) for each chunk. 
if Nchunks > 1
    chunkStart_sec = floor(linspace(args.signalStart,args.signalStart+args.signalLength_sec, Nchunks+1)); %[sec]
    %chunkLength_sec = diff(chunkStart_sec) + [1*ones(1,Nchunks-1), 0]; %[sec]
    chunkLength_sec = diff(chunkStart_sec); %[sec]
    chunkStart_sec = chunkStart_sec(1,1:end-1); %[sec]
else
    chunkStart_sec = args.signalStart; %[sec]
    chunkLength_sec = args.signalLength_sec; %[sec]           
end %if Nchunks > 1
% 
%Convert from sec to samples.
chunkStart_samples = chunkStart_sec * floor(args.fs) + 1; %[samples]
chunkLength_samples = chunkLength_sec * floor(args.fs); %[samples]
%---

%Initialize the variables.
%
%Variables concatenated across the chunks.
allAnomalyPosition  = cell([args.Nchannels, NfreqBands, Nthresholds]);
allEventPosition    = cell([args.Nchannels, Nthresholds]);
allEventLabel       = cell([args.Nchannels, Nthresholds]);
allEventSignal      = cell([args.Nchannels, Nthresholds]);
allEventBPFsignal   = cell([args.Nchannels, Nthresholds]);
allEventT           = cell([args.Nchannels, Nthresholds]);
%
%Memory pre-allocation for speed up the loop.
%
for chk=1:+1:length(chunkStart_samples) %Loop over the chunks.
        
    %Display the number of the chunk being processed.
    disp(['Detection on chunk ' num2str(chk) '/' num2str(length(chunkStart_samples))]);
    
    %Compute the current chunk for all the channels.
    chunk = args.signal(:, chunkStart_samples(chk)+(0:+1:chunkLength_samples(chk)-1));

    %---
    %Reflect the time series.
    %(in order to minimize edge artifacts due to the transient response of the filters).
    %
    %Compute the index for the fraction of the signal to be reflected.
    args.indFraction = round(args.reflectFraction * size(chunk,dimSamples));
    %    
    chunk_r = [chunk(:,args.indFraction:-1:1), chunk, chunk(:,end:-1:end-args.indFraction+1)];
    % %DEBUGGING: Show the signals.
    % figure, plot([nan(size(chunk_r,dimChannels),args.indFraction),...
    %               chunk_r(:,args.indFraction+1:end-args.indFraction),...
    %               nan(size(chunk_r,dimChannels),args.indFraction)].','.-b')
    % hold on, plot(chunk_r.','o-r')
    %
    %Ref: 
    %Mike X. Cohen, Analyzing Neural Time Series Data, Theory and Practice, MIT press, 2014.
    %Figure 7.3, p 78
    %14.9 Filtering Each Trial versus Filtering Concatenated Trials, p 190 
    %Figure 14.11, p 191

    %Compute the settling time of the filters.
    %Valid if and only if the signal was reflected.
    args.indSettling = args.indFraction + 1;                
    %---
    
    %Reset the variables from the previous chunk and initialize the dimensions.
    %NOTE: This initialization ensures that the anomalyPos variable has the
    %proper dimensions even in the case some frequency bands were not processed 
    %because they fail to satisfy the Nyquist criterion.
    anomalyPos(1:args.Nchannels,1:NfreqBands,1:Nthresholds) = {[]};
    %
    %Memory pre-allocation for speed up the loop.
    chunk_bpf = NaN(args.Nchannels,size(chunk,dimSamples),NfreqBands);
    %ampThreshold_high = NaN(args.Nchannels,NfreqBands,Nthresholds);
    %ampThreshold_low = NaN(args.Nchannels,NfreqBands,Nthresholds);
    ampThreshold = NaN(args.Nchannels,NfreqBands,Nthresholds);
    %
    parfor bb=1:NfreqBands %Loop over the frequency bands (parfor).

        %---
        %Extract the parameters for the LFDR configuration (compatible with the parfor).
        %NOTE: The entire structure 'args' is broadcasted to all the
        %workers. Be aware of the possible communication overhead.
        LFDR = args.LFDR;
        %---

        %---
        %Check if the sampling rate satisfy the Nyquist criterion.
        if args.fs <= (2*args.BPFcfg{bb}.f2)
            %Continue with the next frequency band.
            continue,
        end %if args.fs <= (2*args.BPFcfg{bb}.f2)
        %---

        %---
        %Band-pass Filtering.
        FDFout = function_FDF_v1(chunk_r.', args.BPFcfg{bb});
        %---
        
        %---
        %Restore the length of the signals to remove the transient of the filtering.
        FDFout.filteredSignal = FDFout.filteredSignal.';
        chunk_bpf(:,:,bb) = FDFout.filteredSignal(:,args.indSettling:end-(args.indSettling-1));
        %---
        
        %---
        %Normalize the band-pass filtered signals: Spectral whitening.
        %Z-score normalization of each channel (rows) across the samples (columns).
        chunk_bpf(:,:,bb) = function_zscore_v1(squeeze(chunk_bpf(:,:,bb)).').';
        %---
        
%         %---
%         %DEBUGGING: Show the resulting signals.
%         cc = 1;
%         figure, plot(squeeze(chunk_bpf(cc,:,bb)),'-b')
%         %---             
        
        %Define an auxiliary variables for compatibility with the parfor.
        %
        bpfSignal = chunk_bpf(:,:,bb);
        %
        %Reset the variables from the previous frequency band (bb-1).
        %
        %ampThreshold_high_ = NaN(args.Nchannels,Nthresholds);
        %ampThreshold_low_ = NaN(args.Nchannels,Nthresholds);        
        ampThreshold_ = NaN(args.Nchannels,Nthresholds);
        %
        anomalyPos_   = {};
        %
        for cc=1:+1:args.Nchannels %Loop over the channels.
        %Refs:
        %https://fr.mathworks.com/help/parallel-computing/parfor.html 
        %https://fr.mathworks.com/help/parallel-computing/nested-parfor-loops-and-for-loops.html
        %IMPORTANT:
        %If you want to reduce parallel overhead and speed up your computation, run the outer loop in parallel.
        %If you convert the inner loop instead, then each iteration of the outer loop initiates a separate parfor-loop.
        %That is, the inner loop conversion creates 100 parfor-loops. Each of the multiple parfor executions incurs overhead.
        %If you want to reduce parallel overhead, you should run the outer loop in parallel instead, because overhead only occurs once.
                    
            %Compute the threshold via LFDR.
            for th=1:+1:Nthresholds %Loop over the LFDR threshold values.
            
            %---
            %Set the LFDR threshold value.
            LFDR.lfdr_threshold = args.LFDR.lfdr_threshold(th);    
                
            %Compute the LFDR.
            %LFDRout = function_LFDR_v1(squeeze(chunk_bpf(cc,:,bb)), LFDR);
            LFDRout = function_LFDR_v1(bpfSignal(cc,:), LFDR);
            
            %ampThreshold_high_(cc,th) = LFDRout.empThreshold_high;
            %ampThreshold_low_(cc,th) = LFDRout.empThreshold_low;
            ampThreshold_(cc,th) = max(abs([LFDRout.empThreshold_high, LFDRout.empThreshold_low]));

            %Refs:
            %https://www.tonmeister.ca/wordpress/2016/12/15/probability-and-death/
            %https://www.tonmeister.ca/wordpress/2017/04/27/probability-density-functions-part-3/
            %https://www.tonmeister.ca/wordpress/2017/04/26/probability-density-functions-part-2/
            %So, what we can see in Figures 7 through 12 (inclusive) is that,
            %regardless of the original PDF of the signal, if you band-limit it, the result has a Gaussian distribution.
            %And yes, I tried other bandwidths and filter slopes. The result, generally speaking, is the same.
            %One part of this effect is a little obvious. The high-pass filter (in this case, at 200 Hz)
            %removes the DC component, which makes all of the PDF’s symmetrical around the 0 line.
            %However, the "punch line" is that, regardless of the distribution of the signal coming into 
            %your system (and that can be quite different from song to song as I showed in this posting) 
            %the PDF of the signal after band-limiting (say, being sent to your loudspeaker drivers) will be Gaussian-ish.
            %---
                
            %---
            %Detection of anomalies (i.e. amplitude outliers) using the LFDR thresholds.
            
            if args.ampEnv(bb)
                %Compute the amplitude envelope.
                %ampEnv = abs(hilbert(squeeze(chunk_bpf(cc,:,bb))));
                ampEnv = abs(hilbert(bpfSignal(cc,:)));
            else
                %Do not compute the amplitude envelope.
                %ampEnv = squeeze(chunk_bpf(cc,:,bb));
                ampEnv = bpfSignal(cc,:);
            end %if args.ampEnv(bb)

            %Compute the indices using the threshold on the amplitude envelope.
            indTh = ampEnv <= ampThreshold_(cc,th);
            %
            %DEBUGGING: Show the resulting signals.
            %figure, plot(ampEnv)
            %
            %Apply the threshold on the amplitude envelope.
            ampEnv(indTh) = 0;
            %
            %DEBUGGING: Show the resulting signals.
            %hold on, plot(ampEnv)

            %Find the coordinates of the anomalies (peaks after thresholding the amplitude envelope).
            [~, anomaliesInd] = findpeaks(ampEnv,'MinPeakWidth',0,'MaxPeakWidth',Inf,'MinPeakDistance',(args.windowLength/2)*args.fs);

            %Refs:
            %https://fr.mathworks.com/help/signal/ug/peak-analysis.html
            %https://fr.mathworks.com/help/signal/ref/findpeaks.html
            %---
            
            %---
            if isempty(anomaliesInd)
                anomalyPos_{cc,th} = [];
            else
                %Compute the anomaly time position in seconds.
                anomalyPos_{cc,th} =...
                (chunkStart_samples(chk) + anomaliesInd - 1) / args.fs; %[sec]
            end %if isempty(anomaliesInd)
            %---
            
            end %Loop over the LFDR threshold values.    
            
        end %Loop over the channels.

        %---
        %Update the variables in a way compatible with the parfor.
        %
        %ampThreshold_high(:,bb,:) = ampThreshold_high_;
        %ampThreshold_low(:,bb,:) = ampThreshold_low_;            
        ampThreshold(:,bb,:) = ampThreshold_; 
        %
        anomalyPos(:,bb,:) = anomalyPos_;
        %---

%         %DEBUGGING: Show the resulting signals (band-pass filtered and denoised),
%         %the thresholds and the markers for the anomalies.
%         %
%         cc = 1;
%         %
%         if args.ampEnv(bb)
%             %ampEnv = abs(hilbert(squeeze(chunk_bpf(cc,:,bb))));
%             ampEnv = abs(hilbert(bpfSignal(cc,:)));
%         else
%             %ampEnv = squeeze(chunk_bpf(cc,:,bb));
%             ampEnv = bpfSignal(cc,:);
%         end %if args.ampEnv(bb)
%         %
%         figure, hold on,
%         plot(chunk(cc,:),'-k','LineWidth',1)
%         %plot(squeeze(chunk_bpf(cc,:,bb)),'-b','LineWidth',2)
%         plot(bpfSignal(cc,:),'-b','LineWidth',2)
%         plot(ampEnv,':k','LineWidth',1)
%         plot(anomaliesInd,ampEnv(anomaliesInd),'+r','MarkerSize',20)
%         %
%         axis tight, haxes = gca;
%         %
%         %plot([haxes.XLim(1), haxes.XLim(end)],[ampThreshold_high(cc,bb,th), ampThreshold_high(cc,bb,th)],'--k','LineWidth',1)
%         %plot([haxes.XLim(1), haxes.XLim(end)],[ampThreshold_low(cc,bb,th), ampThreshold_low(cc,bb,th)],'--k','LineWidth',1)        
%         plot([haxes.XLim(1), haxes.XLim(end)],[ampThreshold(cc,bb,th), ampThreshold(cc,bb,th)],'--k','LineWidth',1)
%         %  
%         %plot(abs(chunk(cc,:))+ampThreshold(cc,bb,th),'--k','LineWidth',1)
        
     end %Loop over the frequency bands (parfor).
     
     %---
     %Concatenate the markers structure fields across the chunks.
     if chk == 1 %If che == 1, then initialize.
         allAnomalyPosition = anomalyPos;
     else %If chk > 1, then concatenate.
         allAnomalyPosition{:,:,:} = [allAnomalyPosition{:,:,:}, anomalyPos{:,:,:}];
     end %if chk == 1
     %---      
     
%      %---
%      %DEBUGGING: Show the resulting signals.
%      cc = 1;
%      figure, hold on,
%      plot(squeeze(chunk_bpf(cc,:,1)),'-b')
%      plot(squeeze(chunk_bpf(cc,:,2)),'-k')
%      plot(squeeze(chunk_bpf(cc,:,3)),'-r')   
%      %---
     
     %---
     %Merge the anomalies co-ocurrent in time across the frequency bands, to obtain a single event.
     %
     %IMPORTANT: Note that here we merge anomalies which are co-occurrent in time
     %but pertain to different frequency bands. Her we do not merge anomalies co-occurring
     %in time within the same frequency bands. As a consequence, the resulting events
     %can still be superimposed in time (i.e. occuring within the same finite length time window).
     %If further co-occurrence constraints are required, co-occurent events should be identified 
     %and merged in a single event outside this function.

     %Reset the variables from the previous chunk.
     eventPos   = {};
     eventLabel = {};
     %
     sortedEventPos   = {};
     sortedEventLabel = {};
     %
     %Memory pre-allocation for speed up the loop.
     %
     for th=1:+1:Nthresholds %Loop over the LFDR threshold values.
         
     for cc=1:+1:args.Nchannels %Loop over the channels.
         for bb1=1:+1:NfreqBands %First loop over the frequency bands.
             
             %Initialize the time position of the events with the time position
             %of the anomalies in the current frequency band (bb1).
             eventPos{cc,bb1,th}   = anomalyPos{cc,bb1,th};
             %
             %Initialize the variables.
             eventLabel{cc,bb1,th} = [];
             for aa=1:+1:length(anomalyPos{cc,bb1,th}) %Loop over the anomalies.
                 
                 %---
                 %Clear the label from the previous anomaly (aa-1).
                 label = zeros(1,NfreqBands);
                 %
                 %Initialize the label of the event with the label of the 
                 %anomaly aa in the frequency band bb1.
                 label(bb1) = 1; %Logical label for the anomaly aa in the frequency band bb1.
                 for bb2=bb1+1:+1:NfreqBands %Second loop over the frequency bands.

                     %Determine if any anomaly in the frequency band bb2 is within
                     %the time window of the current anomaly pertaining to the frequency band bb1.
                     indStartAnomaly = anomalyPos{cc,bb2,th} > (anomalyPos{cc,bb1,th}(aa)-(args.windowLength/2));
                     indEndAnomaly = anomalyPos{cc,bb2,th} < (anomalyPos{cc,bb1,th}(aa)+(args.windowLength/2));
                     %
                     indAnomaly = indStartAnomaly & indEndAnomaly;
                     
                     %Update the parameters of the current event.
                     if sum(indAnomaly) >= 1
                         %
                         %Update the position of the event (i.e. merged anomalies).
                         eventPos{cc,bb1,th}(aa) = mean([anomalyPos{cc,bb1,th}(aa), anomalyPos{cc,bb2,th}(indAnomaly)]);
                         %
                         %Update the label of the event (i.e. merged anomalies).
                         label(bb2) = 1; %Logical label corresponding the co-ocurrent anomaly in the frequency band bb2.
                         %
                         %Clear the anomaly(ies) of the frequency band bb2 co-occurring 
                         %with the anomaly of the frequency band bb1, since they
                         %were merged in a single event. This prevents anomaly(ies)
                         %processed in the inner loop bb2 from being considered again
                         %in the outer loop bb1.
                         anomalyPos{cc,bb2,th}(indAnomaly) = [];
                         %
                     end %if sum(indAnomaly) >= 1                     

                 end %Second loop over the frequency bands.                     
                 %---

                 %---
                 %Compute the label characterizing the current event (i.e. merged anomalies).
                 eventLabel{cc,bb1,th}{aa} = strrep(num2str(label),' ','');
                 %--- 
                     
             end %Loop over the anomalies.
         end %First loop over the frequency bands.

         %Remove the frequency band index, since now it is meaningless due
         %to the fact the anomalies have been merged across the frequency bands.
         catEventPos   = cat(dimCols,[eventPos{cc,:,th}]);
         catEventLabel = cat(dimCols,[eventLabel{cc,:,th}]);
         %
         %Sort the events according to ascending position.
         [sortedEventPos{cc,th}, indSort] = sort(catEventPos,'ascend');
         sortedEventLabel{cc,th} = catEventLabel(indSort);
         %
         %Refs:
         %https://fr.mathworks.com/help/matlab/ref/double.cat.html
         %https://fr.mathworks.com/help/matlab/ref/sort.html
         
         %Concatenate the markers structure fields across the chunks.
         allEventPosition{cc,th}   = [allEventPosition{cc,th},   sortedEventPos{cc,th}];
         allEventLabel{cc,th}      = [allEventLabel{cc,th},      sortedEventLabel{cc,th}];
         
     end %Loop over the channels.
     
     end %Loop over the LFDR threshold values.
     %---   

     %---
     %Use the saveSignal and saveBPFsignal flag to record the time series
     %corresponding to each detected event, only if they were asked as 
     %outputs in the current call to this function.
     if saveSignal || saveBPFsignal     
     
     for cc=1:args.Nchannels %Loop over the channels.
     %Refs:
     %https://fr.mathworks.com/help/parallel-computing/parfor.html 
     %https://fr.mathworks.com/help/parallel-computing/nested-parfor-loops-and-for-loops.html
     %IMPORTANT:
     %If you want to reduce parallel overhead and speed up your computation, run the outer loop in parallel.
     %If you convert the inner loop instead, then each iteration of the outer loop initiates a separate parfor-loop.
     %That is, the inner loop conversion creates 100 parfor-loops. Each of the multiple parfor executions incurs overhead.
     %If you want to reduce parallel overhead, you should run the outer loop in parallel instead, because overhead only occurs once.        
             
         parfor th=1:Nthresholds %Loop over the LFDR threshold values (parfor).
     
         %Compute the number of samples in the time window.
         halfWindowLength_samples = round(args.windowLength * args.fs / 2);
         %
         for ee=1:+1:length(sortedEventPos{cc,th}) %Loop over the events.
             %Compute the indices for the time window.
             indt = (round(sortedEventPos{cc,th}(ee)*args.fs)-halfWindowLength_samples:+1:round(sortedEventPos{cc,th}(ee)*args.fs)+halfWindowLength_samples);
                   
             %Check for valid temporal indices.
             if indt(1) >= 1 && indt(end) <= size(chunk,dimSamples)
                 
                 %Time vector for the time series around each detected event.
                 allEventT{cc,th} = [allEventT{cc,th}, {(indt + chunkStart_samples(chk) - 1) / args.fs}]; %[sec]
             
                 %Extract the raw time series around each detected event.
                 if saveSignal
                     allEventSignal{cc,th} = [allEventSignal{cc,th}, {chunk(cc,indt)}];
                 end %if saveSignal
             
                 %Extract the band-pass time series around each detected event.
                 if saveBPFsignal
                     allEventBPFsignal{cc,th} = [allEventBPFsignal{cc,th}, {squeeze(chunk_bpf(cc,indt,:))}];
                 end %if saveBPFsignal               
                 
             else
                       
                 %Time vector for the time series around each detected event.
                 allEventT{cc,th} = [allEventT{cc,th}, {[]}]; %[sec]
             
                 %Extract the raw time series around each detected event.
                 if saveSignal
                     allEventSignal{cc,th} = [allEventSignal{cc,th}, {[]}];
                 end %if saveSignal
             
                 %Extract the band-pass time series around each detected event.
                 if saveBPFsignal
                     allEventBPFsignal{cc,th} = [allEventBPFsignal{cc,th}, {[]}];
                 end %if saveBPFsignal
             
             end %if indt(1) >= 1 && indt(end) <= size(chunk,dimSamples)
               
         end %Loop over the events.

         end %Loop over the LFDR threshold values (parfor).

     end %Loop over the channels.
     
     end %if saveSignal || saveBPFsignal
     %---

end %Loop over the chunks.
%--------------------------------------------------------------------------

%Synthesize the output structure ------------------------------------------

%Remove the signal field.
if isfield(args,'signal')
    args = rmfield(args,'signal');
end %if isfield(args,'signal')

%Synthesize the markers structure.
markers = struct(...
'anomalyPosition', {allAnomalyPosition},...
'eventPosition',   {allEventPosition},...
'eventLabel',      {allEventLabel});

Dout = struct(...
'markers',      markers,...
'ampThreshold', ampThreshold,...
'Dcfg',         args);    

Sout = struct(...
'eventSignal', {allEventSignal},...
'eventT',      {allEventT});

BPFout = struct(...
'eventBPFsignal', {allEventBPFsignal},...
'eventT',         {allEventT},...
'eventF',         {args.BPFcfg});

%--------------------------------------------------------------------------

end %Main function.

%% checkInputs function.

function args = checkInputs(signal, varargin)

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
%validSignal = @(x) isnumeric(x) && size(x,2)>size(x,1) && ~any(any(isnan(x)));
validSignal = validArray;
p.addRequired('signal', validSignal);

%Add optional name-value pair argument into input parser scheme.
validChannelLabel = @(x) (iscellstr(x));
defaultChannelLabel = {}; %Consecutive channel numbers will be assigned after parsing (see below).
p.addParameter('channelLabel', defaultChannelLabel, validChannelLabel);
%
validReflectFraction = @(x) isnumeric(x) && isscalar(x) && ~isnan(x) && x >= 0;
defaultReflectFraction = 0; %Fraction of the signal to be reflected before filtering.
p.addParameter('reflectFraction', defaultReflectFraction, validReflectFraction);
%
validFs = validScalar;
defaultFs = []; %This is non a valid sampling rate value (an error message will be generated below).
p.addParameter('fs', defaultFs, validFs);
%
validSignalStart = validScalar;
defaultSignalStart = 0; %[sec]
p.addParameter('signalStart', defaultSignalStart, validSignalStart);
%
validWindowLength = validScalar;
defaultWindowLength = 0.2; %[sec]
p.addParameter('windowLength', defaultWindowLength, validWindowLength);
%
validAmpEnv = @(x) islogical(x);
%Flags to compute the amplitude envelope in each frequency band.
defaultAmpEnv = false(size(varargin{1,1}{:}.BPFcfg));
p.addParameter('ampEnv', defaultAmpEnv, validAmpEnv);
%
validBPFcfg = @(x) (isstruct(x) || iscell(x));
defaultBPFcfg = {struct()}; %The fields of this structure are checked inside the filtering function (function_FDF_v1.m).
p.addParameter('BPFcfg', defaultBPFcfg, validBPFcfg);
%
validLFDR = @(x) (isstruct(x));
defaultLFDR = struct(); %The fields of this structure are checked below and inside the LFDR function (function_LFDR_v1.m).
p.addParameter('LFDR', defaultLFDR, validLFDR);
%
%NOTE: We desided not to use the positional argument format.
%For instance:
%
% %Add optional, positional argument into input parser scheme.
% validFs = validScalar;
% defaultFs = [];
% p.addOptional('fs', defaultFs, validFs);

%ADD MORE INPUT PARAMETERS HERE...
%---

%---
%Parse the input arguments.
%This code line supports the input arguments be included in a single 
%structure (i.e. Dcfg), or in the name-value pair format.
p.parse(signal, varargin{1,1}{:});

%Return the valid arguments.
args = p.Results;
%---

%---
%Define specific constrains on the input arguments.

dimRows = 1; %Rows.
dimCols = 2; %Columns.

%Compute the default labels for the channels.
if isempty(args.channelLabel)
    args.channelLabel = cellstr(strsplit(num2str( (1:+1:size(args.signal,dimRows)) )));
end %isempty(args.channelLabel)

%Check the compatibility between the channels and their labels.
assert(size(args.signal,dimRows) == length(args.channelLabel),...
       ['The number of channels of the input signal and the number of labels are not the same.']);
   
%Check for a valid sampling rate value.
assert(~isempty(args.fs),...
       ['A non valid sampling rate value was detected.']);

%Check for the size of the ampEnv array.   
assert(length(args.ampEnv)==length(args.BPFcfg),...
       ['Non valid size for the ampEnv array.']);

%Check for a valid LFDR threshold values.   
if isfield(args.LFDR,'lfdr_threshold') || ~isempty(args.LFDR.lfdr_threshold)
    assert(~any(isnan(args.LFDR.lfdr_threshold(:))),...
           ['Non valid LFDR threshold values detected.']);
else
    %Default value: Allow 10 percent of false positive outliers of amplitude.
    %args.LFDR.lfdr_threshold = 0.1;
    
    error('MATLAB:function_node',...
          ['LFDR threshold values not defined.']);

    %Refs: 
    %\Papers\statistical analysis\FalseDiscoveryRate\Efron2004.pdf
end %if isfield(args.LFDR,'lfdr_threshold') || ~isempty(args.LFDR.lfdr_threshold)
   
if args.reflectFraction == 0
    warning('MATLAB:function_node',...
            ['Possible filtering transient artifacts can be produced since the input signal is not reflected before filtering.']);
end %if args.reflectFraction == 0

%ADD MORE SPECIFIC CONSTRAINS HERE...
%---

%Refs:
%https://www.mathworks.com/help/matlab/ref/nargin.html
%https://www.mathworks.com/help/matlab/ref/varargin.html
%https://fr.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html
%https://www.mathworks.com/help/matlab/ref/inputparser.html            
            
end %checkInputs function.

%% Departamento de Fisica Medica, Instituto Balseiro, Centro Atomico Bariloche, Argentina.
%https://fisica.cab.cnea.gov.ar/fisicamedica/sites-ddellavale/

%Authors: Damian Dellavale (dellavale@cab.cnea.gov.ar, dellavaledamian@gmail.com).
%Date: 18/06/2018 (updated 15/10/2020).

%% Description.

%In this script the Frequency Domain Filtering (FDF) is implemented.
%If f0 and Bw or f1 and f2 are specified, a Band-Pass filter is implemented.
%If just f2 (or f1=NaN) is specified, a Low-Pass filter is implemented.
%If just f1 (or f2=NaN) is specified, a High-Pass filter is implemented.
%
%The filter design procedure is described in:
%/matlab_functions/filtering/Frequency_Domain_Filtering/misNotas/FDFnotes.docx

%% Tree of dependencies.

%function_FDF_v1.m
% function_window_v1.m

%% Changes from previous versions.

%None.

%NOTE:
%Premature optimization is the root of all evil ! 
%Donald Knuth (https://en.wikiquote.org/wiki/Donald_Knuth)

%% References.

%\damian\softwareTools_compilers_IDEs\MATLAB\matlab_functions\filtering\Frequency_Domain_Filtering

%\damian\tesis-materias-escuelas_IB\tesis Maestria Valencia Urbina
%\damian\projects\imaging\photoplethysmography\2015-2016\scripts\finales22d11d2016

% Examples for "windowParam" configuration:
%
% XdB = -30; %Attenuation at the gaussian window extreme values (tails).
% windowParam = struct('name','gausswin','XdB',XdB,'alpha',sqrt(-2*XdB/(10*log10(exp(1))))); -> Power.
% windowParam = struct('name','gausswin','XdB',XdB,'alpha',sqrt(-XdB/(10*log10(exp(1)))));   -> Amplitude.
%
% windowParam = struct('name','hann','sflag','symmetric');
%
% windowParam = struct('name','hamming','sflag','symmetric');
%
% windowParam = struct('name','rec'); %Rectangular window.
%
%Ref: See Parameters characterizing the window functions at the end of function_window_v1.m

%% Main function.

function [FDFout] = function_FDF_v1(signal, FDFcfg)
%==========================================================================
%Inputs:
%signal -> Data [signal] (matrix: samples x channels).
%          If "signal=NaN", the "indSettling" and "FilterMag" are computed,
%          but the filter is not applied to the signal.
%FDFcfg -> Frequency Domain Filtering configuration (structure array).
%          'f0'   -> Center frequency of the BPF [Hz] (scalar).
%          'Bw'   -> Bandwidth of the Band-Pass Filter [Hz] (scalar).
%          'f1'   -> Lower cutoff frequency (null amplitude = -Inf dB) [Hz] (scalar).
%          'f2'   -> Higher cutoff frequency (null amplitude = -Inf dB) [Hz] (scalar).
%          'type' -> Define the filter as one of the following:
%                    type='stopband', a stoppass filter is implemented.
%                    otherwise (type='bandpass'), a bandpass filter is implemented.
%          'zeropadding' -> Padding flag [none] (scalar).
%                    If zeropadding>=0, the input signal is zero-padded to the "(zeropadding+1)-th" next power of 2 greater than Ns (length of "signal").
%                    If zeropadding<0, the input signal is truncated at the "zeropadding-th" previous power of 2 lesser than Ns (length of "signal").
%          'conv' -> Convolution flag (string).
%                    If conv='linear', zero-padding is implemented so the product of the FFTs results to the LINEAR convolution in the time domain.
%                    Otherwise, no zero-padding is implemented so the product of the FFTs results to the CIRCULAR convolution in the time domain. 
%          'causal' -> Causal filtering flag (logical).
%                    If causal=true, Causal filtering is implemented. That is, the filter kernel (h) is shifted to the rigth so h=0 for all t<0.
%                    Otherwise, Non-causal filtering. The filter kernel (h) is centered at t=0.
%                    IMPORTANT: The 'causal' flag applies only in the case of linear convolution (conv='linear').
%          'freqWindowParam' -> Parameters of the window function in the frequency domain (structure array).
%          'timeWindowParam' -> Parameters of the window function in the time domain (structure array).
%                    'name'  -> Name defining the window type (string).
%                    'alpha' -> Parameter for the gausswin: alpha is inversely proportional to the standard
%                               deviation of a Gaussian random variable.
%                    'sflag' -> Sampling parameter for hamming and hann windows (string).
%                               'sflag' can be either 'periodic' or 'symmetric' (the default).
%                               The 'periodic' flag is useful for DFT/FFT purposes, such as in spectral analysis.
%                               The DFT/FFT contains an implicit periodic extension and the periodic flag enables a signal windowed
%                               with a periodic window to have perfect periodic extension. When 'periodic' is specified,
%                               hamming/hann computes a length "windowLength+1" window and returns the first "windowLength" points.
%                               When using windows for filter design, the 'symmetric' flag should be used. 
%          'Nf'   -> Number of frequencies to evaluate the Filter's frequency response.
%          'compFreqResp' -> Flag to compute the frequency response of the filter (logical).
%          'compSettTime' -> Flag to compute the settling time of the filter (logical).
%          'fs' -> Sampling rate [Hz] (scalar).

%Outputs:
%FDFout -> Output structure (structure).
%          'filteredSignal', Filtered signal [signal] (matrix: samples x channels). 
%          'indSettling', Index corresponding to the output settling time (scalar).
%          'FilterMag', Magnitude of the Filter's frequency response (array). 
%          'f', Frequency vector corresponding to the BPF's frequency response (array).
%          'hlength', Relative number of samples in which h fall under 1 percent of its maximum value (scalar).
%==========================================================================

%Argument completion ------------------------------------------------------
if (nargin < 2) || isempty(signal) || isempty(FDFcfg),
    error('MATLAB:function_FDF','Input argument error".');
end
%--------------------------------------------------------------------------

%Check the input arguments ------------------------------------------------
% if ~isnan(signal(1,1)) && size(signal,1) <= size(signal,2),
%    warning('MATLAB:function_FDF',...
%            'The input signals must be column arrays of the input matrix.') 
% end
%
if ~isnumeric(signal) || any(isnan(signal(:))),
   warning('MATLAB:function_FDF',...
           'The input signals must be a valid matrix of column arrays.') 
end

assert(isstruct(FDFcfg), 'Input argument error in function "function_FDF": FDFcfg must be a structure array.');

if ~isfield(FDFcfg, 'fs') || isempty(FDFcfg.fs) || ~isnumeric(FDFcfg.fs) || ~isfinite(FDFcfg.fs),
    error('MATLAB:function_FDF','Error in the filter configuration (FDFcfg).');
end

if isfield(FDFcfg, 'f0') && isfield(FDFcfg, 'Bw'),
    if length(FDFcfg.f0)==1 && length(FDFcfg.Bw)==1,
        %Compute the cutoff frequencies.
        f1 = FDFcfg.f0 - FDFcfg.Bw/2;
        f2 = FDFcfg.f0 + FDFcfg.Bw/2;
        FilterType = 'BPF'; %BandPass Filter.
    else
        error('MATLAB:function_FDF','Error in the filter configuration (FDFcfg).');
    end
elseif isfield(FDFcfg, 'f1') && isnumeric(FDFcfg.f1) && isfinite(FDFcfg.f1) &&...
       isfield(FDFcfg, 'f2') && isnumeric(FDFcfg.f2) && isfinite(FDFcfg.f2),
    if length(FDFcfg.f1)==1 && length(FDFcfg.f2)==1,
        %Compute the cutoff frequencies.
        f1 = FDFcfg.f1;
        f2 = FDFcfg.f2;
        FilterType = 'BPF'; %BandPass Filter.
    else
        error('MATLAB:function_FDF','Error in the filter configuration (FDFcfg).');        
    end
elseif isfield(FDFcfg, 'f1') && isnumeric(FDFcfg.f1) && isfinite(FDFcfg.f1),
    if length(FDFcfg.f1)==1,
        %Compute the cutoff frequencies.
        f1 = FDFcfg.f1;
        f2 = FDFcfg.fs/2; %[Hz] Nyquist frequency.
        FilterType = 'HPF'; %HighPass Filter.
    else
        error('MATLAB:function_FDF','Error in the filter configuration (FDFcfg).');        
    end    
elseif isfield(FDFcfg, 'f2') && isnumeric(FDFcfg.f2) && isfinite(FDFcfg.f2),
    if length(FDFcfg.f2)==1,
        %Compute the cutoff frequencies.
        f1 = 0; %[Hz]
        f2 = FDFcfg.f2;
        FilterType = 'LPF'; %LowPass Filter.
    else
        error('MATLAB:function_FDF','Error in the filter configuration (FDFcfg).');        
    end    
else
    error('MATLAB:function_FDF','Error in the filter configuration (FDFcfg).');
end

if ~isfield(FDFcfg, 'compFreqResp') || ~islogical(FDFcfg.compFreqResp),
    FDFcfg.compFreqResp = false; %Do not compute the frequency response of the filter.
end
if FDFcfg.compFreqResp && (~isfield(FDFcfg, 'Nf') || isempty(FDFcfg.Nf) || ~isnumeric(FDFcfg.Nf) || ~isfinite(FDFcfg.Nf)),
    FDFcfg.Nf = 2^10; %Default value for the number of frequencies to evaluate the BPF's frequency response.
end

if ~isfield(FDFcfg, 'compSettTime') || ~islogical(FDFcfg.compSettTime),
    FDFcfg.compSettTime = false; %Do not compute the settling time of the filter.
end

if ~isfield(FDFcfg, 'conv') || ~ischar(FDFcfg.conv),
    FDFcfg.conv = 'circular'; %Default value.
    warning('MATLAB:function_FDF','Circular convolution is implemented by default.') 
end

if ~isfield(FDFcfg, 'causal') || ~islogical(FDFcfg.causal),
    FDFcfg.causal = false; %Default value.
    warning('MATLAB:function_FDF','Non-causal convolution is implemented by default.') 
end

if f1<0 || f2>FDFcfg.fs/2 || f1>=f2,
    error('MATLAB:function_FDF','Error in the filter configuration (FDFcfg).');
end

if ~isfield(FDFcfg, 'type') || ~ischar(FDFcfg.type),
    FDFcfg.type = 'bandpass'; %A BPF is implemented.
end
%--------------------------------------------------------------------------

%Default values of the outputs --------------------------------------------
filteredSignal = NaN;
indSettling = NaN;
FilterMag = NaN;
fmag = NaN;
hlength = NaN;

f1aux = f1;
f2aux = f2;
%--------------------------------------------------------------------------

%Parameters ---------------------------------------------------------------
%The input signals are column arrays, so we compute the fft over the Rows.
dim = 1; %Rows.

Nch = size(signal,3-dim); %Number of channels.
Ns = size(signal,dim); %Number of samples.

%Length of the FFT:
%If zeropadding>=0, is the "(zeropadding+1)-th" next power of 2 greater than Ns (length of "signal").
%If zeropadding<0, is the "zeropadding-th" previous power of 2 lesser than Ns (length of "signal").
nfft = 2^(nextpow2(Ns) + FDFcfg.zeropadding); %[samples]
%nfft = 2^(ceil(log(Ns)/log(2)) + FDFcfg.zeropadding); %[samples]

%The length of the onsided PSD depend on the nfft:
%onesidedLength = (nfft/2 + 1) for nfft even.
%onesidedLength = (nfft + 1)/2 for nfft odd -> (never happens). 
onesidedLength = nfft/2 + 1 - mod(nfft,2)/2; %[samples]

%Change the length of the signal to the corresponding power of 2,
%in order to use the fft() matlab function.
if FDFcfg.zeropadding < 0, %If nfft < Ns. 
    signal = signal(1:nfft,:); %Truncation.
else %If nfft >= Ns.
    signal = [signal; zeros(nfft-Ns, Nch)]; %Zero-padding.
end

%Compute the parameters for customized windows.
windowName = FDFcfg.freqWindowParam.name;
switch lower(FDFcfg.freqWindowParam.name), %Switch for window's selection.

    case 'gausswin', %Gaussian window:
        
        FDFcfg.freqWindowParam.alpha = sqrt(-FDFcfg.freqWindowParam.XdB/(10*log10(exp(1))));
        %FDFcfg.freqWindowParam.name = 'gausswin';
        
    case 'gausswinwide1', %Wide Gaussian window:
        
        FDFcfg.freqWindowParam.alpha = sqrt(-FDFcfg.freqWindowParam.XdB/(10*log10(exp(1))));
        f1aux = ( f1 + f2 - (f2-f1)*FDFcfg.freqWindowParam.alpha/sqrt(log(2)) ) / 2; %Cutoff frequency at -3dB.
        f2aux = ( f1 + f2 + (f2-f1)*FDFcfg.freqWindowParam.alpha/sqrt(log(2)) ) / 2; %Cutoff frequency at -3dB.
        FDFcfg.freqWindowParam.name = 'gausswin';
        
    case 'gausswinwide2', %Wide Gaussian window:
        
        FDFcfg.freqWindowParam.alpha = sqrt(-FDFcfg.freqWindowParam.XdB/(10*log10(exp(1))));
        f1aux = ( f1 + f2 - (f2-f1)*FDFcfg.freqWindowParam.alpha ) / 2; %Cutoff frequency at the inflection point.
        f2aux = ( f1 + f2 + (f2-f1)*FDFcfg.freqWindowParam.alpha ) / 2; %Cutoff frequency at the inflection point.
        FDFcfg.freqWindowParam.name = 'gausswin';       
        
    case 'hannwide1', %Wide Hann window:
        
        f1aux = f1 - (f2-f1)/(pi-2); %Cutoff frequency at -3dB.
        f2aux = f2 + (f2-f1)/(pi-2); %Cutoff frequency at -3dB.
        FDFcfg.freqWindowParam.name = 'hann';   
        
    case 'hannwide2', %Wide Hann window:
        
        f1aux = (3*f1 - f2)/2; %Cutoff frequency at the inflection point.
        f2aux = (3*f2 - f1)/2; %Cutoff frequency at the inflection point.
        FDFcfg.freqWindowParam.name = 'hann';        
        
    case 'tukeywide', %Wide Tukey window:
        
        f1aux = ( f1*(FDFcfg.freqWindowParam.r-2) + FDFcfg.freqWindowParam.r*f2 ) / (2*(FDFcfg.freqWindowParam.r-1)); %Cutoff frequency at 0dB.
        f2aux = ( f2*(FDFcfg.freqWindowParam.r-2) + FDFcfg.freqWindowParam.r*f1 ) / (2*(FDFcfg.freqWindowParam.r-1)); %Cutoff frequency at 0dB.
        FDFcfg.freqWindowParam.name = 'tukey';        

    otherwise,
        %A standard window is implemented (see function_window).
end %Switch for window's selection.

%Check the values of f1aux and f2aux.
switch lower(FilterType), %Switch for filter's type.
    
    case 'bpf', %Band-Pass Filter:
        
        if f1aux>=0 && f2aux<=FDFcfg.fs/2,
            f1 = f1aux;
            f2 = f2aux;
        else
            disp(['Frequencies out of range for ',windowName,' window.',char(10),...
                  'f1aux=',num2str(f1aux),'; f2aux=',num2str(f2aux),char(10),...
                  FDFcfg.freqWindowParam.name,' implemented instead with:',char(10),...
                  'f1=',num2str(f1),'; f2=',num2str(f2),char(10)])
        end        
        
    case 'hpf', %High-Pass Filter:
        
        if f1aux>=0,
            f1 = f1aux;
        else
            disp(['Frequency out of range for ',windowName,' window.',char(10),...
                  'f1aux=',num2str(f1aux),char(10),...
                  FDFcfg.freqWindowParam.name,' implemented instead with:',char(10),...
                  'f1=',num2str(f1),char(10)])
        end        

    case 'lpf', %Low-Pass Filter:
        
        if f2aux<=FDFcfg.fs/2,
            f2 = f2aux;
        else
            disp(['Frequency out of range for ',windowName,' window.',char(10),...
                  'f2aux=',num2str(f2aux),char(10),...
                  FDFcfg.freqWindowParam.name,' implemented instead with:',char(10),...
                  'f2=',num2str(f2),char(10)])
        end        

    otherwise,
        %None.
end %Switch for filter's type.   
%--------------------------------------------------------------------------

%Compute the Frequency response of the filter -----------------------------
if FDFcfg.compFreqResp,

FDFmag = FDFcfg; %Make a copy to avoid changes in the original structure.   
    
fmag = linspace(f1,f2,FDFmag.Nf); %Compute the frequency vector.

%Compute the window function.
switch lower(FilterType), %Switch for filter's type.
    case 'bpf', %Band-Pass Filter:
        
        if strcmpi(FDFmag.type,'stopband'),
            FilterMag = 1 - function_window_v1(FDFmag.freqWindowParam, FDFmag.Nf);    
        else
            FilterMag = function_window_v1(FDFmag.freqWindowParam, FDFmag.Nf);
        end
        
    case 'hpf', %High-Pass Filter:
        
        %Scaling the window prameters.
        if isfield(FDFmag.freqWindowParam, 'alpha'),
            FDFmag.freqWindowParam.alpha = FDFmag.freqWindowParam.alpha / 2; %gausswin
        end
        if isfield(FDFmag.freqWindowParam, 'r'),
            FDFmag.freqWindowParam.r = FDFmag.freqWindowParam.r / 2; %tukey
        end
        %Compute the window doubling the number of samples.
        FilterMag = function_window_v1(FDFmag.freqWindowParam, 2*FDFmag.Nf);
        FilterMag = FilterMag(1:length(FilterMag)/2);
        
    case 'lpf', %Low-Pass Filter:
        
        %Scaling the window prameters.
        if isfield(FDFmag.freqWindowParam, 'alpha'),
            FDFmag.freqWindowParam.alpha = FDFmag.freqWindowParam.alpha / 2; %gausswin
        end
        if isfield(FDFmag.freqWindowParam, 'r'),
            FDFmag.freqWindowParam.r = FDFmag.freqWindowParam.r / 2; %tukey
        end
        %Compute the window doubling the number of samples.
        FilterMag = function_window_v1(FDFmag.freqWindowParam, 2*FDFmag.Nf);
        FilterMag = FilterMag(length(FilterMag)/2+1:end);
        
    otherwise,
        
        %None.
        
end %Switch for filter's type.

% %DEBUGGING: Plot the signals.
% figure, plot(fmag, FilterMag,'-b')

end  %if FDFcfg.compFreqResp,
%--------------------------------------------------------------------------

%Compute the transient response of the filter -----------------------------
if FDFcfg.compSettTime,

if abs( mean(signal,1) ./ max(abs(signal),[],1) ) > 1e-6*ones(1,size(signal,2)),
    error('MATLAB:function_FDF',...
          'The input signals must have zero mean, otherwise the settling time is assumed to be inf.')
else
    indSettling = round(10*Ns/100); %10 percent of the samples.
end

end %if FDFcfg.compSettTime,
%--------------------------------------------------------------------------

if ~isnan(signal(1,1)), %If the "signal" is a single NaN, do not apply the filter.

%Compute the frequency vector ---------------------------------------------
f = linspace(0,+FDFcfg.fs/2,onesidedLength).'; f = [f(1:end); -f(end-1:-1:2)];
%--------------------------------------------------------------------------

%Compute the window function ----------------------------------------------

%Compute the indices for locating the window function.
indf1 = find(f>=f1); indf1 = indf1(1); 
indf2 = find(f>=f2); indf2 = indf2(1);

%In case of f2=FDFcfg.fs/2, the following is necessary because the negative frequencies
%(second half of the fft vector) do not include the FDFcfg.fs/2 frequency.
indf2Neg = find(f<=-f2);
if f2 >= FDFcfg.fs/2 || isempty(indf2Neg),
    indf2Neg = onesidedLength;
else
    indf2Neg = indf2Neg(end);
end

%Compute the window length.
windowLength = length(f(indf1:indf2));

%Compute the window function.
switch lower(FilterType), %Switch for filter's type.
    case 'bpf', %Band-Pass Filter:
        
        windowFunction = function_window_v1(FDFcfg.freqWindowParam, windowLength);
        
    case 'hpf', %High-Pass Filter:
        
        %Scaling the window prameters.
        if isfield(FDFcfg.freqWindowParam, 'alpha'),
            FDFcfg.freqWindowParam.alpha = FDFcfg.freqWindowParam.alpha / 2; %gausswin
        end
        if isfield(FDFcfg.freqWindowParam, 'r'),
            FDFcfg.freqWindowParam.r = FDFcfg.freqWindowParam.r / 2; %tukey
        end
        %Compute the window doubling the number of samples.
        windowFunction = function_window_v1(FDFcfg.freqWindowParam, 2*windowLength);
        windowFunction = windowFunction(1:length(windowFunction)/2);
        
    case 'lpf', %Low-Pass Filter:
        
        %Scaling the window prameters.
        if isfield(FDFcfg.freqWindowParam, 'alpha'),
            FDFcfg.freqWindowParam.alpha = FDFcfg.freqWindowParam.alpha / 2; %gausswin
        end
        if isfield(FDFcfg.freqWindowParam, 'r'),
            FDFcfg.freqWindowParam.r = FDFcfg.freqWindowParam.r / 2; %tukey
        end
        %Compute the window doubling the number of samples.
        windowFunction = function_window_v1(FDFcfg.freqWindowParam, 2*windowLength);
        windowFunction = windowFunction(length(windowFunction)/2+1:end);
        
    otherwise,
        
        %None.
        
end %Switch for filter's type.

%Compute the filter in frequency domain.
if strcmpi(FDFcfg.type,'stopband'),
    windowFunction = 1 - windowFunction;
    H = ones(size(signal,dim),1);
else
    H = zeros(size(signal,dim),1);
end
H(indf1:indf1+windowLength-1) = windowFunction;
H(indf2Neg:indf2Neg+windowLength-1) = flipud(windowFunction);

%In case of f1=0Hz or f2=FDFcfg.fs/2, this is necessary because the negative frequencies
%(second half of the fft vector) do not include the 0Hz and FDFcfg.fs/2 frequencies.
H = H(1:nfft);

% %DEBUGGING: Plot the signals.
% figure, plot(f, H,'.b')
%--------------------------------------------------------------------------

%Windowing in time domain (Window method for FIR filter design) -----------

%Compute the impulse response of the filter.
h = ifft(H,nfft,dim); %The zero-time component is in the center of the array.
%
%Check for the residual imaginary part.
if any(abs(imag(h(:))) > 1e-9) && isreal(signal)
    error('MATLAB:function_FDF','The impulse response of the filter has non null imaginary part.');
else
    h = real(h);
end %if any(abs(imag(h(:))) > 1e-9) && isreal(signal)

%Re-arrange h by shifting the zero-time component to the left of the array. 
h = fftshift(h); 
%
%Apply the window in time domain.
win = function_window_v1(FDFcfg.timeWindowParam, nfft);
h = h.*win;
    
%Re-arrange h by shifting the zero-time component back to the center of the array.
h = fftshift(h); 

%---
%If the h values goes to zero in a small numer of samples relative to the
%length of the signal, then we can assume "circular convolution = linear convolution"

%Compute the maximum of h and the corresponding index.
[hmax, ind_hmax] = max(h);

%Compute the indices for the values of h below the 1 percent of its maximum value.
hlength = find( h(ind_hmax:end) < 1*hmax/100 );

%Compute the number of samples between the maximum of h and the first
%sample under 1 percent of its maximum value.
hlength = hlength(1) + ind_hmax - 1;

%Compute the relative number of samples in which h fall under 1 percent of its maximum value.
hlength = hlength / length(h);
%---

%--------------------------------------------------------------------------

%Switch for linear convolution in time domain -----------------------------
if strcmp(FDFcfg.conv,'linear'),
   
    %Update the signal's length required so the product of the FFTs results to the LINEAR convolution in time domain.
    nfft = 2*nfft; %Minimum required length = size(signal,dim) + size(h,dim) - 1 = 2*nfft - 1;
    onesidedLength = nfft/2 + 1 - mod(nfft,2)/2; %[samples]
    
    %Update the frequency vector.
    f = linspace(0,+FDFcfg.fs/2,onesidedLength).'; f = [f(1:end); -f(end-1:-1:2)];
    
    %Apply the zero-padding on the filter kernel.
    if FDFcfg.causal, %Causal filter.
        
        %Re-arrange h by moving the zero-time component to the left of the array. 
        h = fftshift(h); 
        
        h(nfft) = 0; %Zero-pad h to make its length equals to nfft.
        
    else %Non-causal filter.

        Nh = length(h); %Compute the kernel length.
        
        h1 = h(1:Nh/2); %Extract the first half of the kernel.
        h2 = h(Nh/2+1:end); %Extract the last half of the kernel.
        
        %Apply the zero-padding to the first half of the kernel.
        h1(nfft/2) = 0;
        
        %Apply the zero-padding to the last half of the kernel.
        h2 = flipud(h2); 
        h2(nfft/2) = 0;
        h2 = flipud(h2); 
        
        %Reconstruct the zero-padded filter kernel.
        h = [h1; h2];
        
    end
    
    %Apply the zero-padding on the input signal.        
    signal(nfft,:) = 0; %Zero-pad signal to make its length equals to nfft.
   
%     %DEBUGGING: Plot the signals.
%     figure, plot(signal,'.k')
%     figure, plot(h,'.b')
    
end %if strcmp(FDFcfg.conv,'linear'),
%--------------------------------------------------------------------------

%Compute the fft of the filter impulse response h -------------------------
H = fft(h,nfft,dim); 

%Normalize the frequency response.
H = H / max(abs(H)); 

% %DEBUGGING: Plot the signals.
% figure, plot(f, abs(H),'.:b')
% figure, plot(f, angle(H),'.:r')
%--------------------------------------------------------------------------

%Compute the fft of the signal --------------------------------------------
FFT = fft(signal,nfft,dim);
%--------------------------------------------------------------------------

%Reshape "H" to get dimensions to match those of "signal" -----------------
%Thus, in matrix "HH" the column "H" is replicated over the channels ("Nch" times).

HH = H(:,ones(1,Nch)); 

%--------------------------------------------------------------------------

%Filtering in frequency domain --------------------------------------------

%Apply the window in the frequency domain.
filteredFFT = FFT .* HH;

%--------------------------------------------------------------------------

% %DEBUGGING: Plot the signals ----------------------------------------------
% figure, plot(f,abs(FFT),'.-b')
% hold on, plot(f, abs(filteredFFT),'o-r')
% legend('|FFT|','|filteredFFT|','location','northeast')
% legend('boxon')
% 
% box on, grid off,
% axis tight;
% %axis([]);
% 
% set(gca,'xscale','linear');
% set(gca,'yscale','linear');
% %set(gca,'xscale','log');
% %set(gca,'yscale','log');
% 
% set(gca,'fontname','Times New Roman','fontsize',16','fontWeight','bold');
% set(gca,'TickDirMode','manual','TickDir','in');
% set(gca,'YMinorTick','On','XMinorTick','On','Linewidth',1);
% 
% xlabel(['Frequency [Hz]'],'interpreter','LaTex')
% ylabel(['Amplitude [arb. units]'],'interpreter','LaTex')
% %--------------------------------------------------------------------------

%Recover the filtered signal ----------------------------------------------
filteredSignal = ifft(filteredFFT,nfft,dim);
%
%Check for the residual imaginary part.
if any(abs(imag(filteredSignal(:))) > 1e-9) && isreal(signal),
    error('MATLAB:function_FDF','The filtered signal has non null imaginary part.');
else
    filteredSignal = real(filteredSignal);
end %if any(abs(imag(filteredSignal(:))) > 1e-9) && isreal(signal),
%--------------------------------------------------------------------------

%In the case of non-causal linear filtering, recover the original signal
%length by removing the zero-padding --------------------------------------

%Frequency response of nfft samples. 
%H1 = fft(filteredSignal);

if ~FDFcfg.causal && (nfft > Ns),
    filteredSignal = filteredSignal(1:Ns,:);
end %if ~FDFcfg.causal && (nfft > Ns),
%We have verified that the first Ns = length(signal) samples corresponds to
%the non-causal linear convolution. See the linear convolution with a Kronecker's delta
%in the script "conv_circular_vs_linear.m".
%Refs:
%/matlab_functions/filtering/Frequency_Domain_Filtering/misNotas/
%FDFnotes.docx
%conv_circular_vs_linear.m

%Frequency response of Ns samples. 
%H2 = fft(filteredSignal);

%IMPORTANT:
%We have eliminated the discontinuity in "filteredSignal" introduced by the
%zero-padding before the computation of H2. As a consequence,
%H2 show less oscillations in the pass-band (Gibbs effect) that H1.
%--------------------------------------------------------------------------

end %if ~isnan(signal(1,1)), %If the "signal" is a single NaN, do not apply the filter.

%Synthesize the output structure ------------------------------------------

FDFout = struct(...
         'filteredSignal', filteredSignal,...
         'indSettling', indSettling,...
         'FilterMag', FilterMag,...
         'f', fmag,...
         'hlength', hlength);

%--------------------------------------------------------------------------


end %Main function.

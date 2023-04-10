%% Departamento de Fisica Medica, Instituto Balseiro, Centro Atomico Bariloche, Argentina.

%Authors: Damian Dellavale (dellavale@cab.cnea.gov.ar, dellavaledamian@gmail.com).
%Date: 01/01/2016

%% Description.

%In this script we synthesize the weighting function for windowing in time domain.

%% Tree of dependencies.

%None.

%% Changes from previous versions.

%None.

%NOTE:
%Premature optimization is the root of all evil ! 
%Donald Knuth (https://en.wikiquote.org/wiki/Donald_Knuth)

%% References.

%Matlab -> Help -> window
%http://www.mathworks.com/help/signal/ref/window.html
%https://en.wikipedia.org/wiki/Window_function

%% Main function.

function windowFunction = function_window_v1(windowParam, windowLength)
%==========================================================================
%Inputs:
%windowParam    -> Parameters of the window function (structure array).
%                  'name'  -> Name defining the window type (string).
%                  'alpha' -> Parameter for the gausswin: alpha is inversely proportional to the standard
%                             deviation of a Gaussian random variable.
%                  'sflag' -> Sampling parameter for hamming and hann windows (string).
%                             'sflag' can be either 'periodic' or 'symmetric' (the default).
%                             The 'periodic' flag is useful for DFT/FFT purposes, such as in spectral analysis.
%                             The DFT/FFT contains an implicit periodic extension and the periodic flag enables a signal windowed
%                             with a periodic window to have perfect periodic extension. When 'periodic' is specified,
%                             hamming/hann computes a length "windowLength+1" window and returns the first "windowLength" points.
%                             When using windows for filter design, the 'symmetric' flag should be used. 
%                  'r'     -> A Tukey window is a rectangular window with the first and last 100*r/2 percent of the samples equal to parts of a cosine.
%                             r is a real number between 0 and 1. If you input r = 0, you obtain a rectwin window. If you input r = 1, you obtain a hann window.
%windowLength   -> Length of the window function.
%
%Outputs:
%windowFunction -> Synthesized window function.
%==========================================================================

%Argument completion ------------------------------------------------------
if (nargin < 2)||isempty(windowParam)||isempty(windowLength),...
   error('MATLAB:function_window','Input argument error.');
end
%--------------------------------------------------------------------------

%Check the input arguments ------------------------------------------------
assert(isstruct(windowParam), 'Input argument error in function "function_window": windowParam must be a structure array.');
%--------------------------------------------------------------------------

switch lower(windowParam.name), %Switch for window's selection.
    case 'gausswin', %Gaussian window:
        
        windowFunction = gausswin(windowLength, windowParam.alpha);
        %windowFunction = window(@gausswin, windowLength, windowParam.alpha);
        
%         warning('MATLAB:function_window', ['The Gaussian window will modify the DC component of the signal',...
%                                            '(use a rectangular window to leave the DC level unchanged).']);

    case 'hamming', %Hamming window:
        
        windowFunction = hamming(windowLength, windowParam.sflag);
        %windowFunction = window(@hamming, windowLength, windowParam.sflag);
        
%         warning('MATLAB:function_window', ['The Hamming window will modify the DC component of the signal',...
%                                            '(use a rectangular window to leave the DC level unchanged).']);
                                       
    case 'hann', %Hann (Hanning) window:
        
        windowFunction = hann(windowLength, windowParam.sflag);
        %windowFunction = window(@hann, windowLength, windowParam.sflag);
        
%         warning('MATLAB:function_window', ['The Hann window will modify the DC component of the signal',...
%                                            '(use a rectangular window to leave the DC level unchanged).']);   

    case 'tukey', %Tukey (tapered cosine) window:
        
        windowFunction = tukeywin(windowLength, windowParam.r);
        %windowFunction = window(@tukeywin, windowLength, windowParam.r);
        
%         warning('MATLAB:function_window', ['The Tukey window will modify the DC component of the signal',...
%                                            '(use a rectangular window to leave the DC level unchanged).']); 

    otherwise, %Rectangular window:
        %This function is provided for completeness; a rectangular window is equivalent to no window at all.
        windowFunction = rectwin(windowLength);
        %windowFunction = ones(windowLength,1);
end %Switch for window's selection.

% %DEBUGGING: Show the window function in time and frequency domains.
% wvtool(windowFunction) 
% 
% figure, plot(10*log10(windowFunction/max(windowFunction)))
% box on, grid off, axis tight;
% set(gca,'xscale','linear');
% set(gca,'yscale','linear');
% set(gca,'fontname','Times New Roman','fontsize',16','fontWeight','bold');
% set(gca,'TickDirMode','manual','TickDir','in');
% set(gca,'YMinorTick','On','XMinorTick','On','Linewidth',1);
% xlabel('$\mathbf{\textit{index}}$','interpreter','LaTex')
% ylabel('$\mathbf{\textit{Amplitude}~[dB]}$','interpreter','LaTex')
% title(['Window function: ',lower(windowParam.name)],'interpreter','none')

%Parameters characterizing the window functions (extracted from the wvtool):
%1) Leakage factor: ratio of power in the sidelobes to the total window power.
%2) Relative sidelobe attenuation in frequency domain: difference in height from the mainlobe peak to the highest sidelobe peak.
%3) Mainlobe width (-3dB): width of the mainlobe at 3 dB below the mainlobe peak.
%4) Attenuation at the window extreme values (tails), in the time domain.

%Rectangular window:
%1) 9.28 [%]
%2) -13.3 [dB]
%3) 6.8665e-5 [x pi rad/sample]
%4) 0 [dB]

%Hann (Hanning) window:
%1) 0.05 [%]
%2) -31.5 [dB]
%3) 0.00011444 [x pi rad/sample]
%4) -78 [dB]

%Gaussian window ( alpha = sqrt(-2*XdB/(10*log10(exp(1)))), XdB=-78dB ):
%1) ~0 [%]
%2) -175.7 [dB]
%3) 0.00025177 [x pi rad/sample]
%4) -78 [dB]

%Gaussian window ( alpha = sqrt(-2*XdB/(10*log10(exp(1)))), XdB=-15dB ):
%1) 0.01 [%]
%2) -46.1 [dB]
%3) 0.00010681 [x pi rad/sample] => Similar to Hann window.
%4) -15 [dB]

%Gaussian window ( alpha = sqrt(-2*XdB/(10*log10(exp(1)))), XdB=-30dB ):
%1) ~0 [%]
%2) -78.2 [dB]
%3) 0.00015259 [x pi rad/sample]
%4) -30 [dB]

%Gaussian window ( alpha = sqrt(-2*XdB/(10*log10(exp(1)))), XdB=-10dB ):
%1) 0.1 [%]
%2) -34.2 [dB]
%3) 9.1553e-5 [x pi rad/sample]
%4) -10 [dB]

%Hamming window:
%1) 0.04 [%]
%2) -42.7 [dB]
%3) 9.9182e-5 [x pi rad/sample]
%4) -10.97 [dB]

end %Main function.

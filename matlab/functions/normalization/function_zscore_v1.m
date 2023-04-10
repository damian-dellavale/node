%% Departamento de Fisica Medica, Instituto Balseiro, Centro Atomico Bariloche, Argentina.

%Authors: Damian Dellavale (dellavale@cab.cnea.gov.ar, dellavaledamian@gmail.com).
%Date: 14/12/2015

%% Description.

%In this script the input signal is z-score normalized in order to have
%zero mean and unit variance.
%
%The z-score normalization is applied on each channel (columns) taken separately.

%% Tree of dependencies.

%None.

%% Changes from previous versions.

%None.

%% References.

%Detection of epileptiform activity in EEG signals (Gajic, 2015).pdf
%The temporal structures and functional significance of scale-free brain activity, SI (He, 2010).pdf
%High gamma power is phase-locked to theta oscillations, SOM.pdf
%F:\damian_F\projects\EBSP\DBSproject\results\Mesbah_consultation_2015
%https://en.wikipedia.org/wiki/Standard_score

%https://www.mathworks.com/help/stats/zscore.html

%% Main function.

function [normSignal] = function_zscore_v1(rawSignal)
%==========================================================================
%Inputs:
%rawSignal -> Raw signal to be normalized [signal] (column array: samples x channels).

%Outputs:
%normSignal -> Normalized signal [signal] (column array: samples x channels).
%==========================================================================

%Argument completion ------------------------------------------------------
if (nargin < 1), error('MATLAB:function_zscore','Input argument error.'); end
%--------------------------------------------------------------------------

%Check the input arguments ------------------------------------------------
%assert(max(size(rawSignal))==size(rawSignal,1), 'Input argument error in function "function_zscore": The signal must be a column array.');

if (max(size(rawSignal))>size(rawSignal,1)) 
    warning('MATLAB:function_zscore',...
            ['The normalization is performed across the rows. '...
             'The input matrix has Nrows < Ncol. '...
             'Consider traspose the input matrix to implement the normalization across the columns.']);
end %if (max(size(rawSignal))>size(rawSignal,1))

if any(isnan(rawSignal(:)))
    warning('MATLAB:function_zscore',...
            ['The NaN''s detected in the input signal are omited for z-score normalization.']);
end %if any(isnan(rawSignal(:)))
%--------------------------------------------------------------------------

%Default values of the outputs --------------------------------------------
%--------------------------------------------------------------------------

%Parameters ---------------------------------------------------------------
%Set the dimention along to which compute the z-score normalization.
dim = 1; %Rows -> Samples.
%--------------------------------------------------------------------------

%Normalization in order to have zero mean and unit variance ---------------

%Zero mean normalization.
offset = mean(rawSignal,dim,'omitnan');
%offset = repmat(offset,size(rawSignal,dim),1);
offset = offset(ones(size(rawSignal,dim),1),:); %This one is faster, try -> 
%F:\damian_F\softwareTools_compilers_IDEs\MATLAB\matlab_templates\bsxfun_repmat_kron.m
normSignal = rawSignal - offset;

%Unit variance normalization.
%Compute the standard deviation.
%Standard deviation normalized with "N-1", where "N" is the number of samples. 
%See eq. (5) in "Introduction to the Theory of error, by Yardley Beers". 
flag = 0;
stdDev  = std(normSignal,flag,dim,'omitnan');
%stdDev = repmat(stdDev,size(normSignal,dim),1);
stdDev = stdDev(ones(size(normSignal,dim),1),:); %This one is faster, try -> 
%F:\damian_F\softwareTools_compilers_IDEs\MATLAB\matlab_templates\bsxfun_repmat_kron.m
normSignal = normSignal./stdDev;

%Verify the normalizations:
%mean(normSignal,dim)
%var(normSignal,dim)
%--------------------------------------------------------------------------

end %Main function.


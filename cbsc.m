function EstimatedSourceCount = cbsc( MicSigs, AnalysisParams )
%CBSC Estimate the number of sources from two microphone signals
%
% Syntax:	EstimatedSourceCount = ...
%                          cbsc( MicSigs, AnalysisParams )
%
% Inputs:
%          MicSigs - The two microphone signals at a node. Samples are in
%                    each row, microphones in each column and nodes in each
%                    page.
% 	AnalysisParams - The structure of analysis parameters:
%   |-> .Nfft            - FFT samples.
%   |-> .ovlap           - Percentage of overlap.
%   |-> .segLen          - Segment length in samples.
%   |-> .MaxClusters     - Maximum number of sources to count.
%   |-> .ClustDist       - Distance Metric.
%   |-> .ClustEval       - Evaluation criterion.
%   '-> .ClustReplicates - Number of times to repeat clustering using new
%                          initial cluster centroid positions.
%
% Outputs:
% 	EstimatedSourceCount - The estimated number of sources.
%
% Example:
% 	x_2ch = [];               % Microphone signals
% 	SrcCount = cbsc( x_2ch ); % Count
%
% See also: mscohere, kmeans, evalclusters

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 1 September 2017
% Version: 1.0 (1 September 2017)
% Version: 0.1 (21 April 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
  AnalysisParams.Nfft = 256;                     % FFT samples
  AnalysisParams.ovlap = 0.5;                    % Overlap percentage
  AnalysisParams.segLen = 32000;                 % Segment length samples
  AnalysisParams.MaxClusters = 6;                % Maximum source count
  AnalysisParams.ClustDist = 'correlation';      % Distance Metric
  AnalysisParams.ClustReplicates = 100;          % Number of repeats
  AnalysisParams.ClustEval = 'CalinskiHarabasz'; % Evaluation criterion 
end

%%
Nnodes = size(MicSigs,3);

EstimatedSourceCount=[];
for node = 1:Nnodes   
    
    
    %% Compute Magnitude-Squared Coherence 
    win = hamming(AnalysisParams.Nfft,'p');
    noverlap = AnalysisParams.Nfft * AnalysisParams.ovlap;
    
    Nseg = AnalysisParams.segLen;
    msc=[];
    for k = 1:(size(MicSigs,1)/Nseg)
                
            m1seg = MicSigs( (k-1)*Nseg+1:k*Nseg, 1, node);
            m2seg = MicSigs( (k-1)*Nseg+1:k*Nseg, 2, node);
            
            msc(k,:) = mscohere( ...
            m1seg, m2seg, ...
                win, noverlap, AnalysisParams.Nfft );

    end    
    
    %% Cluster MSC values
    clust = zeros( size(msc,1), AnalysisParams.MaxClusters );
    for i = 1:AnalysisParams.MaxClusters
        clust(:,i) = kmeans( msc, i, ...
            'Distance', AnalysisParams.ClustDist, ...
            'Replicates', AnalysisParams.ClustReplicates );
    end
    
    %% Evaluate clusters and estimate source count
    eva = evalclusters( msc, clust, AnalysisParams.ClustEval);
    EstimatedSourceCount(node) = eva.OptimalK;
    
end
end

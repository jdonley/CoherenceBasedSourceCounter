
%% Settings
runInParallel = true;                          % Set true runs processes in parallel, false runs in serial

Nsources = 4;                                  % Number of sources
Nnodes = 1;                                    % Number of nodes

L = [4 7 3];                                   % Room dimensions [x y z] (m)
T = [1 3 0];                                   % Table Size [x y z] (m)
S = T + [0.5 0.5 0]*2;                         % Seating Area [x y z] (m)
Tloc = L/2;                                    % Table Location [x y z] (m)
InterSrcDist  = 0.5;                           % Minimum distance between talkers (metres)
IntraNodeDist = 0.1;                           % Distance between microhpones (metres)

% Set synthesis parameters
SynthParams.SpeechDir = 'Speech_Files\';       % Directory to read speech samples from
SynthParams.c = 343;                           % Sound velocity (m/s)
SynthParams.res = 100;                         % Samples per metre
SynthParams.SNR = 40;                          % AWGN at <SNR>dB
SynthParams.ReverbTime = 0.2;                  % seconds
SynthParams.RIRlength = 1.0;                   % seconds
SynthParams.ForceNodeLoc = Tloc;               % Forces node location (comment out for randomised location)

% Set analysis parameters
AnalysisParams.Nfft = 256;                     % samples
AnalysisParams.ovlap = 0.5;                    % percentage
AnalysisParams.segLen = 32000;                 % samples
AnalysisParams.MaxClusters = 6;                % Maximum number of sources to count
AnalysisParams.ClustDist = 'correlation';      % Distance Metric
AnalysisParams.ClustReplicates = 100;          % Number of times to repeat clustering using new initial cluster centroid positions
AnalysisParams.ClustEval = 'CalinskiHarabasz'; % Evaluation criterion





%% Loop over a variable to determine its influence on Success Rate (SR)
% Initialise (maximum of two variables)
variables = {'Nsources', 'SynthParams.ReverbTime' };              % Variables to test
varnames  = { 'Sources', 'Reverberation Time (RT60) (seconds)' }; % Names of variables being tested
vals      = {       2:6,              0.2:0.1:0.8 };              % Values of the variables

% variables = {'Nsources', 'SynthParams.SNR' };                     % Variables to test
% varnames  = { 'Sources', 'Signal to Noise Ratio (SNR) (dB)' };    % Names of variables being tested
% vals      = {       2:6,            0:5:40 };                     % Values of the variables

% variables = {'Nsources', 'IntraNodeDist' };                       % Variables to test
% varnames  = { 'Sources', 'Intra-node Distance (metres)' };        % Names of variables being tested
% vals      = {       2:6,       0.1:0.1:1 };                       % Values of the variables

Niters    = 100;                                                  % Number of iterations to test each variable



%% Loop
Successes=zeros(Niters,numel(vals{1}),numel(vals{2})); % Initialise matrix
parforArg = int64(runInParallel * realmax);            % Number of workers
if runInParallel, ppool=gcp; end                       % Create parpool
fprintf('\t Completion: '); n=0; tic;                  % Start loop timer

for valInd1 = 1:numel(vals{1})
    val1 = vals{1}(valInd1);
    eval([variables{1} '=val1;']);
    
    for valInd2 = 1:numel(vals{2})
        val2 = vals{2}(valInd2);
        eval([variables{2} '=val2;']);
        
        parfor (iter = 1:Niters, parforArg)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Perform Random Recording
            [x_2ch, fs] = getReceivedSignals( ...
                Nsources, Nnodes, ...
                L, T, S, Tloc, ...
                InterSrcDist, IntraNodeDist, ...
                SynthParams );
            
            %% Count the number of sources
            SrcCount = cbsc( x_2ch, AnalysisParams );
            
            %% Save Successes
            Successes(iter,valInd1,valInd2) = Nsources==SrcCount;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end
        n = showTimeToCompletion( ((valInd1-1)*numel(vals{2}) + valInd2)/(numel(vals{1})*numel(vals{2})), n);
    end
end
delete(ppool); % Delete Parallel Pool

SR  = squeeze(sum(Successes,1) / Niters * 100).'; % Compute Success Rate (SR)


%% Plot Success Rate (SR) results

markLst = {'+','o','.','^','x','s','d','*','v','>','<','p','h'};
lineLst = {'- ','--',': ','-.'};

methodtxt = {'MSC'};

Nmethods = numel(methodtxt);
Nsrc = numel(vals{1});

fH = figure(1);
fH.Name = ['Success Rate - ' variables{1} '_vs_' variables{2}];
ax = gca;
ax.ColorOrder = [0, 0, 0];
set(0,'DefaultAxesColorOrder', [0,0,0]);

ax.LineStyleOrder = ...
    cell2mat([reshape(markLst(repmat(1:Nsrc,Nmethods,1))',[],1), ...
    reshape(lineLst(repmat(1:Nmethods,Nsrc,1)),[],1)]);
set(0,'DefaultAxesLineStyleOrder', cell2mat([reshape(markLst(repmat(1:Nsrc,Nmethods,1))',[],1), ...
    reshape(lineLst(repmat(1:Nmethods,Nsrc,1)),[],1)]));

plot(ax, vals{2}, SR ); hold on

h = zeros(Nsrc+Nmethods, 1);
h(1:Nsrc) = plot( NaN(Nsrc) ); set(h(1:Nsrc),'LineStyle', 'none');
h(Nsrc+1:end) = plot( NaN(Nmethods) ); set(h(Nsrc+1:end),'Marker', 'none');
arrayfun(@(i) set(h(Nsrc+i),'LineStyle',lineLst{i}),1:Nmethods);
srcTxt = [' ' varnames{1}];
legend(h, ...
    [mat2cell([num2str([vals{1}]'), repmat(srcTxt,Nsrc,1)],ones(Nsrc,1),numel(srcTxt)+1); ...
    methodtxt'], ...
    'Location', 'eastoutside');

% legend();
limbuff = 0.05; %5 percent
xrange = [min(vals{2}) max(vals{2})];
yrange = [0 100];
xlabel(varnames{2}); xlim(xrange+[-1 1]*(limbuff*diff(xrange)));
ylabel('Success Rate (%)'); ylim(yrange+[-1 1]*(limbuff*diff(yrange)));
hold off

%% Save results and save plot
fname = [variables{1} '_vs_' variables{2} '_' [methodtxt{:}]];
save([fname '.mat'], 'SR', 'vals', 'variables', 'varnames');
print([fname '.png'],'-dpng','-r600');


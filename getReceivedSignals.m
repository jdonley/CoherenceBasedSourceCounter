function [MicSigs, Fs, SrcLocs, NodeLocs] = ...
    getReceivedSignals( ...
    Nsrcs, Nnodes, ...
    RoomSz, TableSz, SeatSz, TableLoc, ...
    InterSrcDist, IntraNodeDist, ...
    SynthParams )
%GETRECEIVEDSIGNALS Synthesises received signals at microphones in a meeting scenario
%
% Syntax:	[MicSigs, Fs, SrcLocs, NodeLocs] = ...
%                 getReceivedSignals( ...
%                 Nsrcs, Nnodes, ...
%                 RoomSz, TableSz, SeatSz, TableLoc, ...
%                 InterSrcDist, IntraNodeDist, ...
%                 SynthParams )
%
% Inputs:
%           Nsrcs - Number of speech sources (talkers)
%          Nnodes - Number of dual microphone nodes
% 	       RoomSz - Room dimensions [x, y, z] in metres
% 	      TableSz - Size of the meeting table [x, y, z] in metres
% 	       SeatSz - Seating area of the talkers [x, y, z] in metres
% 	     TableLoc - Centre location of the table [x, y, z] in metres
% 	 InterSrcDist - Minimum distance between talkers in metres
%   IntraNodeDist - Distance between microphones of a node in metres
% 	  SynthParams - The structure of synthesis parameters:
%       |-> .SpeechDir    - Directory to read speech samples from.
%       |                   (default: 'Speech_Files\')
%       |-> .c            - Sound velocity (metres/ssecond).
%       |                   (default: 343)
%       |-> .res          - Resolution in samples per metre.
%       |                   (default: 100)
%       |-> .SNR          - Signal to noise ration (for awgn function).
%       |                   (default: 40)
%       |-> .ReverbTime   - Reverberation time in seconds.
%       |                   (default: 0.2)
%       |-> .RIRlength    - RIR length in seconds.
%       |                   (default: 1.0)
%       |-> .ForceNodeLoc - Forces node location (do not define for 
%       |                   randomised location).
%       |                   (default: <randomised location>)
%       '-> .Verbose      - False suppresses output to command window.
%                           (true or false, default: true)
%
% Outputs:
% 	 MicSigs - The synthesised micrphone signals
% 	      Fs - Sampling frequency of MicSigs
% 	 SrcLocs - The locations of the talkers
% 	NodeLocs - The location of the nodes
%
% Example:
%     Nsources = 4;           % Number of sources
%     Nnodes = 1;             % Number of nodes
%     L = [4 7 3];            % Room dimensions [x y z] (m)
%     T = [1 3 0];            % Table Size [x y z] (m)
%     S = T + [0.5 0.5 0]*2;  % Seating Area [x y z] (m)
%     Tloc = L/2;             % Table Location [x y z] (m)
%     InterSrcDist  = 0.5;    % Minimum distance between talkers (metres)
%     IntraNodeDist = 0.1;    % Distance between microhpones (metres)
%     [x_2ch, fs] = getReceivedSignals( ...
%         Nsources, Nnodes, ...
%         L, T, S, Tloc, ...
%         InterSrcDist, IntraNodeDist );
%
% Reference:
%   This function makes use of the Room Impulse Response Generator by
%   Emanuel A.P. Habets: https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator 
%   Note: When SynthParams.Verbose is false the copyright notice is not 
%   printed to the command window. The copyright notice for the included 
%   rir_generator C++ code is:
%  Room Impulse Response Generator (Version 2.1.20141124) by Emanuel Habets
%  Copyright (C) 2003-2014 E.A.P. Habets, The Netherlands.
%
% See also: cbsc, rir_generator, getAllFiles

% Author: Jacob Donley
% University of Wollongong
% Email: jrd089@uowmail.edu.au
% Copyright: Jacob Donley 2017
% Date: 8 September 2017
% Version: 1.1 (8 September 2017)
% Version: 1.0 (1 September 2017)
% Version: 0.1 (21 April 2017)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9
    SynthParams = struct;
end

SynthParams = setSPDefaults(SynthParams);

%% First run attempt rir_generator compilation from current directory
compileRIRGenerator;

%%
rng('shuffle');

Nmics = 2;
res = SynthParams.res;

IntraNodeDists = repmat(IntraNodeDist,Nnodes,1);

NodeAngleRand = rand(Nnodes,1)*2*pi;
AvailableNodePos = zeros(RoomSz(1)*res,RoomSz(2)*res);
AvailableSrcPos  = zeros(RoomSz(1)*res,RoomSz(2)*res);

xvec = round((TableLoc(1)-RoomSz(1)/2)*res+1:(TableLoc(1)+RoomSz(1)/2)*res);
yvec = round((TableLoc(2)-RoomSz(2)/2)*res+1:(TableLoc(2)+RoomSz(2)/2)*res);
Txvec = round((TableLoc(1)-TableSz(1)/2)*res+1:(TableLoc(1)+TableSz(1)/2)*res);
Tyvec = round((TableLoc(2)-TableSz(2)/2)*res+1:(TableLoc(2)+TableSz(2)/2)*res);
Sxvec = round((TableLoc(1)-SeatSz(1)/2)*res+1:(TableLoc(1)+SeatSz(1)/2)*res);
Syvec = round((TableLoc(2)-SeatSz(2)/2)*res+1:(TableLoc(2)+SeatSz(2)/2)*res);

[xx,yy]=meshgrid(xvec/res,yvec/res);
[Txx,Tyy]=meshgrid(Txvec/res,Tyvec/res);

AvailableNodePos( Txvec,  Tyvec) = 1;
AvailableSrcPos ( Sxvec,  Syvec) = 1;
AvailableSrcPos ( Txvec,  Tyvec) = 0;

ChosenTInd = randi(sum(AvailableNodePos(:)),Nnodes,1);
ChosenTXLoc = Txx(ChosenTInd);
ChosenTYLoc = Tyy(ChosenTInd);

ChosenSXLoc=[];ChosenSYLoc=[];
for s = 1:Nsrcs
    % ChosenSInd = randi(sum(AvailableSrcPos(:)), Nsrcs,1);
    ChosenSInd = randi(sum(AvailableSrcPos(:)), 1,1);
    xx_=xx.'.*AvailableSrcPos;xx_=xx_(xx_~=0);
    yy_=yy.'.*AvailableSrcPos;yy_=yy_(yy_~=0);
    ChosenSXLoc(s) = xx_(ChosenSInd);
    ChosenSYLoc(s) = yy_(ChosenSInd);
    msk = (xx.'-ChosenSXLoc(s)).^2 + (yy.'-ChosenSYLoc(s)).^2 > InterSrcDist^2;
    AvailableSrcPos = AvailableSrcPos .* msk;
end

% Source Locations
SrcLocs  = [ChosenSXLoc.', ChosenSYLoc.', repmat(TableLoc(3),Nsrcs,1)];

% Node Locations
NodeLocs = [ChosenTXLoc, ChosenTYLoc, repmat(TableLoc(3),Nnodes,1)];
if isfield(SynthParams, 'ForceNodeLoc') && ~isempty(SynthParams.ForceNodeLoc)
    NodeLocs = SynthParams.ForceNodeLoc;
end

% Channel Locations
[xOffset,yOffset] = pol2cart([NodeAngleRand NodeAngleRand], [-1*IntraNodeDists 1*IntraNodeDists]/2 );
ChannelLocs = NodeLocs(reshape([1:Nnodes; 1:Nnodes],1,[]),:) + ...
    [reshape(xOffset.',[],1),reshape(yOffset.',[],1),zeros(Nnodes*2,1)];



%% Prepare source signals
files = getAllFiles( SynthParams.SpeechDir );
files(randperm(numel(files),numel(files)-Nsrcs))=[];

SrcSigs=[];
SrcTimePos=0;
for file = 1:numel(files)
    [s, Fs] = audioread( files{file} );
    SrcTimePos(file+1) = numel(s)+SrcTimePos(file);
    SrcSigs = [SrcSigs, ...
        circshift([ s.'; ...
        repmat(zeros(1,size(s,1)),numel(files)-1,1);],file-1) ];
end
SrcSigs=SrcSigs.';

%%
slen = size(SrcSigs,1)+(SynthParams.RIRlength * Fs)-1;
MicSigs=zeros(slen,Nmics,Nnodes);

for node=1:Nnodes
    
    for src = 1:Nsrcs
        
        rirARGS = { ...
            SynthParams.c, ...
            Fs, ...
            {ChannelLocs(node*2-1,:), ChannelLocs(node*2,:)}, ...
            SrcLocs(src,:), ...
            RoomSz, ...
            SynthParams.ReverbTime, ...
            SynthParams.RIRlength * Fs};
        
        for mic = 1:Nmics
            if SynthParams.Verbose
                MicSigs(:,mic,node) = MicSigs(:,mic,node) ...
                    + fconv( ...
                    SrcSigs(:,src), ...
                    rir_generator( ...
                    rirARGS{1:2}, rirARGS{3}{mic}, rirARGS{4:end}) );
            else
                evalc( [...
                    'MicSigs(:,mic,node) = MicSigs(:,mic,node)' ...
                    '+ fconv(' ...
                    'SrcSigs(:,src),' ...
                    'rir_generator(' ...
                    'rirARGS{1:2}, rirARGS{3}{mic}, rirARGS{4:end}) )' ...
                    ] );
            end
        end
        
    end
    
    for mic = 1:Nmics
        MicSigs(:,mic,node) = awgn( MicSigs(:,mic,node), SynthParams.SNR, 'measured' );
    end
    
end




end

%--------------------------------------------------------------------------

function SP = setSPDefaults(SP) % Set the default Synthesis Parameters (SP) 
if ~isfield(SP,'SpeechDir')
    SP.SpeechDir = 'Speech_Files\'; % Directory to read speech samples from
end
if ~isfield(SP,'c')
    SP.c = 343;                     % Sound velocity (m/s)
end
if ~isfield(SP,'res')
    SP.res = 100;                   % Samples per metre
end
if ~isfield(SP,'SNR')
    SP.SNR = 40;                    % AWGN at <SNR>dB
end
if ~isfield(SP,'ReverbTime')
    SP.ReverbTime = 0.2;            % seconds
end
if ~isfield(SP,'RIRlength')
    SP.RIRlength = 1.0;             % seconds
end
if ~isfield(SP,'Verbose')
    SP.Verbose = true;              % true or false
end
end

%--------------------------------------------------------------------------

function [y]=fconv(x, h)
%FCONV Fast Parallelised Convolution
%   [y] = FCONV(x, h) convolves x and h in the frequency domain
%
%      x = input vector
%      h = input vector
% 
%      See also CONV
%
% Source URL: https://github.com/JacobD10/SoundZone_Tools
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x(:);
h=h(:);

Ly=size(x,1)+size(h,1)-1;  % 
Ly2=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly

if isa(x, 'gpuArray')
    Ly  = gpuArray(Ly);
    Ly2 = gpuArray(Ly2);
end

X=fft(x, Ly2);             % Fast Fourier transform
H=fft(h, Ly2);	           % Fast Fourier transform
Y=X.*H;        	           %
y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
y=y(1:1:Ly,:);             % Take just the first N elements

end

%--------------------------------------------------------------------------

function compileRIRGenerator

if exist('rir_generator') ~= 3 ... % If MEX file doesn't exist
        && exist('rir_generator.cpp', 'file') == 2 ... % and if C++ source exists
        && ~isempty(mex.getCompilerConfigurations('C++')) % and if C++ compiler exists
    mex -setup C++
    mex rir_generator.cpp
end

end
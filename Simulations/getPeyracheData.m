function spikes = getPeyracheData(shift)

if ~exist('shift','var') || isempty(shift)
    shift = 1;
end
% filepath = '/mnt/cortex-data-40/data-AA/Peyrache/mPFC_Data/'; % folders for each session are in this root directory: change here for other paths
% filepath = '/home/gabriel/Matlab/ssaPaper/Peyrache/mPFC_Data/'; % folders for each session are in this root directory: change here for other paths
filepath = 'C:\Users\Gabma\OneDrive\Dokumente\MATLAB\savedServer\2024\ParisServer\Peyrache\mPFC_Data\';

sessionID = '181020'; % only learning session with nNeurons>50 (nNeurons=55)

spkdata = importdata([filepath sessionID '/' sessionID '_SpikeData.dat']);  % [time-stamp (ms); neuron]
spikes = sortrows(spkdata); spikes(:,1) = spikes(:,1)/1000;

behaviour = importdata([filepath sessionID '/' sessionID '_WakeEpoch.dat'])/1000; % start and end times of Wake epoch (ms/1000 -> s)

% % Optional parameters:
% swsPre = importdata([filepath sessionID '/' sessionID '_SwsPRE.dat'])/1000;
% swsPost = importdata([filepath sessionID '/' sessionID '_SwsPOST.dat'])/1000;
% trialInfo = importdata([filepath sessionID '/' sessionID '_Behavior.dat']); trialInfo(:,1:2) = trialInfo(:,1:2)/1000; % [start_time, stop_time, rule, correct, went_left, light_left] where the rule can be go right (1), go to lit arm (2), go left (3), go to unlit arm (4)
% positions = importdata([filepath sessionID '/' sessionID '_Pos.dat']); positions(:,1) = positions(:,1)/1000;

% spikes = Restrict(spikes, behaviour);

if shift == 2
    spikes = Restrict(spikes, behaviour);
end

if shift == 1
    spikes = Restrict(spikes, behaviour, 'shift','on'); % spikes in the awake session
end
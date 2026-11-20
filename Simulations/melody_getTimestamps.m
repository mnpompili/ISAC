function assemblyActivationTimestamps = melody_getTimestamps(cst, sessionEnd)
% Returns the activation timestamps (s) of the 4 assemblies of the melody
% simulations (whichSimulation = 51:62 in the simulation scripts).
% The melody has 8 notes (108 bpm) and each note activates one assembly:
% notes 1-3 -> assembly 1, note 4 -> assembly 2, notes 5-7 -> assembly 3,
% note 8 -> assembly 4. The melody is repeated at a fixed rate during the
% whole session (default 1200 s = 20 minutes, as in the simulation scripts).
%
% ----cst = whichSimulation - 50: number of melodies per minute----
% whichSimulation | 51  52  53  54  55  56  57  58  59  60  61  62
% melodies/min    | 12  11  10   9   8   7   6   5   4   3   2   1
% repetitions     | 240 220 200 180 160 140 120 100  80  60  40  20
%
% Corrected in November 2026: the previous version repeated the melody every
% 5*cst s during the first 10 minutes only (12/cst melodies per minute).

if nargin < 2
    sessionEnd = 1200; % session length (s) of the melody simulations
end
if ~ismember(cst, 1:12)
    error('cst must be 1 to 12 (whichSimulation 51 to 62)')
end

bpm = 108;
secondsPerQuarter = 60/bpm/4;
notedQuarters = [0 1 2 3 8 9 10 11]; % onsets of the 8 notes
baseTime = notedQuarters * secondsPerQuarter;

melodiesPerMinute = 13 - cst; % cst = 1:12 -> 12:-1:1 melodies per minute
numRepeats = round(melodiesPerMinute*sessionEnd/60); % number of melodies in the session
repeats = (0:numRepeats-1) * 60/melodiesPerMinute; % start time of each melody (s)
timestamps = bsxfun(@plus, repeats', baseTime)+0.03;
assemblyActivationTimestamps = {reshape(timestamps(:,[1 2 3]), [],1); timestamps(:,4); reshape(timestamps(:,[5 6 7]),[],1); timestamps(:,8)};

function [] = savefromCNTdata( cntfile )
% Compatible with the EEG data from
% http://sccn.ucsd.edu/~arno/fam2data/publicly_available_EEG_data.html
%
% Reads a .CNT and its corresponding .DAT files from the dataset and uses
% EEGLAB to save the channel values and the events in a .MAT file.
%
% Takes the base file name of the .CNT file as a parameter.
% Example usage: 
%   savefromCNTdata('fsa1ff01');
% This will create a fsa1ff01.mat file.

% Start EEGLAB, load the .CNT file and save the channel values
eeglab;
EEGcnt = pop_loadcnt( [cntfile '.cnt' ], 'dataformat', 'int16');
EEG.data = EEGcnt.data;

% Read the corresponding .DAT file and save the events/stimuli
[events.event events.resp events.type events.correct events.latency] = ...
    textread([cntfile '.dat'],'%d %d %d %d %d','headerlines',20);
for index = 1:length(events.type)
    tmp = num2str(events.type(index));
    if tmp(3) == '0' 
    	EEGcnt.event(index).type = 'Animal';
    elseif tmp(3) == '5'
        EEGcnt.event(index).type = 'Distractor';
    end
end
EEG.event = EEGcnt.event;

% Create the .MAT file and save the channel values and the events
save([cntfile '.mat'],'EEG');
clear all;

end


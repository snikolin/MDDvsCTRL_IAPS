function data = MDDvsCTRL_events(data)

%% Specify trigger-values and marker channel
MARKER      = strcmp(data.label,'Marker');
NEUTRAL     = 254;
NEGATIVE    = 253;
BLANK       = 255;

%% Find trigger values
i               = [];
trigger         = data.trial{1}(MARKER,:);
trigger(trigger == 0) = 1;
trigger(trigger == 247) = 255;
trigger(1)      = 100;                      % set initial trigger to value 100
trigger(2)      = 101;
trigger(end)    = 100;                      % set final trigger to end value 100   
dif_trigger     = [0,diff(trigger)~=0];     % is 1 when value changes
trigger         = trigger.*dif_trigger;     % set trigger to 0 if it doesn't change
i               = find(trigger);            % indices of triggers

event_data      = [];
event_data(:,1) = trigger(i);                    % Value of the amplifier (252:255)
event_data(:,2) = i;                             % Sample location of the events
event_data(:,3) = i/data.fsample;                % Time of event in seconds
event_data(:,4) = [diff(event_data(:,2)); 0];    % Duration of events in samples
event_data(:,5) = event_data(:,4)/data.fsample;  % Duration of events in seconds

short_trls      = event_data(:,5) < 1e-03;
short_trls(1)   = 0;
event_data(short_trls,:) = [];

%% Categorise events 
% NBACK
stim = find(event_data(:,1) == 254);
for q = 1:length(stim)
    if event_data((stim(q) + 1),1) == 252
        event_data(stim(q),6) = 1; % Hit        
    elseif event_data((stim(q) + 1),1) == 255
        event_data(stim(q),6) = 2; % Miss
    end
end

nostim = find(event_data(:,1) == 253);
for r = 1:length(nostim)
    if event_data((nostim(r) - 1),1) == 255 && event_data((nostim(r) + 1),1) == 252
        event_data(nostim(r),6) = 3; % False alarm (FA)     
    elseif event_data((nostim(r) - 1),1) == 255 && event_data((nostim(r) + 1),1) == 255
        event_data(nostim(r),6) = 4; % Stimulus
    end
end

% RESTING STATE
EEG_triggers = find(event_data(:,5) > 285 & event_data(:,5) < 315 & event_data(:,1) == 254);
event_data(EEG_triggers,6) = 8; % Resting state             

% IAPS
% Between EEG_triggers
BLANK       = 255; 
NEUTRAL     = 254; 
NEGATIVE    = 253; 

% Following an adjustment to Inquisit protocol, the above changed to
NEUTRAL2    = 234;
NEGATIVE2   = 235;

for iaps = EEG_triggers(1)+2:EEG_triggers(2)-1
    if event_data(iaps,1) == BLANK && event_data(iaps,5) > 1.5 && event_data(iaps,5) < 2.5
        event_data(iaps,6) = 5; % Blanks    
    elseif (event_data(iaps,1) == NEUTRAL || event_data(iaps,1) == NEUTRAL2) && event_data(iaps,5) > 4.5 && event_data(iaps,5) < 5.5
        event_data(iaps,6) = 6; % Neutral
    elseif (event_data(iaps,1) == NEGATIVE || event_data(iaps,1) == NEGATIVE2) && event_data(iaps,5) > 4.5 && event_data(iaps,5) < 5.5
        event_data(iaps,6) = 7; % Negative
    end
end 

for nback_breaks = EEG_triggers(2) + 1:length(event_data)
    if event_data(nback_breaks,1) == BLANK && event_data(nback_breaks,5) > 2.1
        event_data(nback_breaks,6) = 9;
    end
end

event_data(event_data(:,6) == 0,:) = [];

%% Create Events list
TRIGGERS = 1:9;
TYPE = {'HIT','MISS','FA','STIM',...
    'BLANK','NEUTRAL','NEGATIVE',...
    'RS-EEG','NBACK_BLOCK'};

event_index = find(event_data(:,6));
for n = 1:length(event_index)
    [~,j] = intersect(TRIGGERS,event_data(event_index(n),6));
    data.event(n).type      = TYPE{j};                       % event type (hit etc...)
    data.event(n).value     = event_data(event_index(n),1);  % trigger value
    data.event(n).sample    = event_data(event_index(n),2);  % sample index
    data.event(n).timestamp = event_data(event_index(n),3);  % time of event
    data.event(n).sampdur   = event_data(event_index(n),4);  % duration (samples)
    data.event(n).duration  = event_data(event_index(n),5);  % duration (seconds)
end

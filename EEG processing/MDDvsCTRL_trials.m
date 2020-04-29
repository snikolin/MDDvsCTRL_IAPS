function trials = MDDvsCTRL_trials(data)
% Segments the data into trials 

%% Create EEG trials
EEG             = find(strcmp({data.event.type},'RS-EEG'));
EEG_LABELS      = {'Eyes-Closed EEG','Eyes-Open EEG'};
chs             = 1:length(data.label);

for s = 1:length(EEG)
    tempdata    = data;
    EEG_start   = tempdata.event(EEG(s)).sample;
    EEG_end     = tempdata.event(EEG(s)).sample + tempdata.event(EEG(s)).sampdur - 1;
    index       = EEG_start:EEG_end;

    % Amount of full 1sec segments that can be extracted from trial data
    total_secs  = floor(length(index)/tempdata.fsample);
    remainder   = length(index) - total_secs * tempdata.fsample;
    
    % Ignore first <1s remainder of data and keep the rest
    A = tempdata.trial{1}(chs,index((remainder + 1):end));
    B = mat2cell(A,length(chs),repmat(tempdata.fsample,total_secs,1));
    tempdata.trial          = B;

    C = tempdata.time{1}(1,index((remainder + 1):end));
    D = mat2cell(C,1,repmat(tempdata.fsample,total_secs,1));
    tempdata.time           = D;
    
    for t = 1:length(tempdata.trial)
        tempevent(t).type       = 'eeg';
        tempevent(t).value      = 254;
        tempevent(t).sample     = index(remainder + 1) + ((t - 1) * tempdata.fsample);
        tempevent(t).timestamp  = tempdata.time{t}(1,1);
        tempevent(t).sampdur    = size(tempdata.trial{t},2);
        tempevent(t).duration   = tempevent(t).sampdur/tempdata.fsample;
    end
    
    tempdata.sampleinfo     = [];
    tempdata.sampleinfo(:,1)= index(remainder + 1):tempdata.fsample:EEG_end;
    tempdata.sampleinfo(:,2)= [(tempdata.sampleinfo(2:end,1) - 1); EEG_end];
    tempdata.event          = tempevent;
    tempdata.block          = EEG_LABELS{s};

    trials(s) = tempdata;    
end

%% Create PRACTICE NBACK trials
% Practice, 
PRE_STIM    = 0.5; % seconds
POST_STIM   = 1.5; % seconds

temptrial   = [];
temptime    = [];
sampleinfo  = [];
tempdata    = [];
tempdata    = data;
tempdata.event = [];

% Extract 2s trials
if EEG(1) > 1
    practice_index  = 1:EEG(1)-1;
    practice_trials = find([data.event(practice_index).sample] - round((data.fsample*PRE_STIM)) > 0);
    for d = 1:length(practice_trials)
        index_start     = data.event(practice_trials(d)).sample - round((data.fsample * PRE_STIM));
        index_end       = data.event(practice_trials(d)).sample + round((data.fsample * POST_STIM)) - 1;
        if index_start >= 1
            temptrial{d}    = data.trial{1}(chs,index_start:index_end);
            temptime{d}     = -PRE_STIM:1/data.fsample:(POST_STIM - (1/data.fsample)); 
            sampleinfo(d,:) = [index_start, index_end];
        end
    end
    
    tempdata.event  = data.event(1:EEG(1)-1);
end

tempdata.trial      = temptrial;
tempdata.time       = temptime;
tempdata.sampleinfo = [];
tempdata.sampleinfo = sampleinfo;
tempdata.block      = 'Practice';

trials(3)   = tempdata;

%% IAPS
PRE_STIM    = 1; % seconds
POST_STIM   = 6; % seconds

temptrial   = [];
temptime    = [];
sampleinfo  = [];
tempdata    = [];
tempdata    = data;

i_iaps = find(ismember({data.event.type},{'NEGATIVE','NEUTRAL'}));

% Extract 7s trials
for e = 1:length(i_iaps)
    index_start     = data.event(i_iaps(e)).sample - round((data.fsample * PRE_STIM));
    index_end       = data.event(i_iaps(e)).sample + round((data.fsample * POST_STIM)) - 1;
    temptrial{e}    = data.trial{1}(chs,index_start:index_end);
    temptime{e}     = -PRE_STIM:1/data.fsample:(POST_STIM - (1/data.fsample)); 
    sampleinfo(e,:) = [index_start, index_end];
end

tempdata.trial      = temptrial;
tempdata.time       = temptime;
tempdata.sampleinfo = [];
tempdata.sampleinfo = sampleinfo;
tempdata.block      = 'IAPS';
tempdata.event      = [];
tempdata.event      = data.event(i_iaps);

trials(4)   = tempdata;

%% N-BACK 0-back, 1-back, 2-back
% Find blocks - last i_nback is the end of the session
nback_blocks = find(ismember({data.event.type},'NBACK_BLOCK'));

PRE_STIM    = 0.5; % seconds
POST_STIM   = 1.5; % seconds
nback       = [1,5;5,9;9,13];
NBACK_LABS  = {'0BACK','1BACK','2BACK'};

row = length(trials);

for seg = 1:size(nback,1)
    
    temptrial   = [];
    temptime    = [];
    sampleinfo  = [];
    tempdata    = [];
    
    indx_events = nback_blocks(nback(seg,1))+1:nback_blocks(nback(seg,2))-1;
    indx_nback  = setdiff(indx_events,nback_blocks);
    
    for f = 1:length(indx_nback)
        index_start     = data.event(indx_nback(f)).sample - round((data.fsample * PRE_STIM));
        index_end       = data.event(indx_nback(f)).sample + round((data.fsample * POST_STIM)) - 1;
        temptrial{f}    = data.trial{1}(chs,index_start:index_end);
        temptime{f}     = -PRE_STIM:1/data.fsample:(POST_STIM - (1/data.fsample)); 
        sampleinfo(f,:) = [index_start, index_end];
    end
    
    tempdata            = data;
    tempdata.trial      = temptrial;
    tempdata.time       = temptime;
    tempdata.sampleinfo = [];
    tempdata.sampleinfo = sampleinfo;
    tempdata.event      = data.event(indx_nback);
    tempdata.block      = NBACK_LABS{seg};
    
    trials(row + seg)   = tempdata;
end

%% Structure data file
% Practice (n-back)
% Eyes-Closed
% IAPS
% Eyes-Open
% N-back

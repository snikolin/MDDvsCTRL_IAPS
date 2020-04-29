%% EEG: MDD VS CONTROLS - ANALYSIS PIPELINE
% Comparison of MDD participants with age- and gender-matched healthy
% controls
%
% Task order
%   2-back practice
%   Eyes closed resting-state
%   IAPS
%   Eyes open resting-state
%   0-back
%   1-back
%   2-back
%
% Pipeline steps
%   Load data
%   Pre-processing
%   Generate events
%   Create trials/epochs
%   Data inspection
%   Automated trial/epoch rejection
%   Visual trial/epoch rejection
%   Summary statistic trail/epoch rejection
%   Channel interpolation
%   ICA
%   ICA component inspection
%   ICA component rejection
%   Re-referencing to common average reference
%   Analyses: e.g. ERP/TFR/Resting-State etc...
%
% Folder organisation
%   0.raw
%   1.preprocess
%   2.events
%   3.trials
%   4.autorej
%   5.visrej
%   6.allrej
%   7.interp
%   8.comps
%   9.postica
%   10.reref
%   11.IAPS

%% General Script Format Per Section
DIR_IN  = 'study_directory';
INFOLD  = 'input_folder'; 
OUTFOLD = 'output_folder'; 
VAR     = 'save_variable';

% File containing filenames for all participant sessions
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    
    % Load data
    load(filename);
    
    % Add data transform function in here
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% Load Data
DIR_IN  = 'study_directory';
INFOLD  = 'TMSi.files';
OUTFOLD = '0.raw'; 
VAR     = 'data';

% File containing filenames for all participant sessions
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    file     = participants(p).name;
    filename = fullfile(DIR_IN,INFOLD,strcat(file,'.Poly5'));
    
    % Read data into fieldtrip format from Poly5
    data     = MDDvsCTRL_readTMSi(filename);
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% Pre-processing
DIR_IN  = 'study_directory';
INFOLD  = '0.raw'; 
OUTFOLD = '1.preprocess'; 
VAR     = 'data';

% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    
    % Load data
    load(filename);
    
    % Line Noise Removal/Filtering/Downsampling  
    data = MDDvsCTRL_preprocess(data);
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% Generate Events
DIR_IN  = 'study_directory';
INFOLD  = '1.preprocess'; 
OUTFOLD = '2.events'; 
VAR     = 'events';

% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    
    % Load data
    load(filename);
    
    % Identify events in the continuous data using triggers
    events = MDDvsCTRL_events(data);
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% Create Trials/Epochs
DIR_IN  = 'study_directory';
INFOLD  = '2.events'; 
OUTFOLD = '3.trials'; 
VAR     = 'trials';

% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    
    % Load data
    load(filename);
    
    % Create trials/epochs from continuous data
    trials = MDDvsCTRL_trials(events);
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% Individual Participant Data Inspection
% ICA can detect repeated noise events (e.g. muscle or blink)
% Look for isolated instances of large artefacts (e.g. person touching EEG
% cap), as these may not be cleaned by ICA
DIR_IN  = 'study_directory';
INFOLD  = '3.trials'; 

% Inspect participants one at a time
PARTICIPANT = 82; 

% Participants
load(fullfile(DIR_IN,'participants.mat'))

% Load data
file        = participants(PARTICIPANT).name;
filename    = fullfile(DIR_IN,INFOLD,file);
load(filename);
    
for a = 1:length(autorej)
    
    % Quick visual inspection of all trials - identify possible BAD trials
    cfg             = [];
    cfg.channel     = 1:64;
    cfg.viewmode    = 'butterfly'; % or 'vertical'
    cfg.ylim        = [-100 100];
    cfg.blocksize   = round(trials(a).time{1}(end) - trials(a).time{1}(1));
    cfg.continuous  = 'no';
    cfg.colorgroups = 'sequential'; 
    cfg.layout      = 'easycapM1.mat';
    
    % Note down potential trials for rejection e.g. in spreadsheet
    ft_databrowser(cfg,trials(a))
end

fprintf(file);

%% Automated Trial/Epoch Rejection
DIR_IN  = 'study_directory';
INFOLD  = '3.trials'; 
OUTFOLD = '4.autorej'; 
VAR     = 'autorej';

% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    
    % Load data
    load(filename);
    
    % Automated trial rejection 
    autorej = MDDvsCTRL_trialrejection(trials);
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% Visual Trial/Epoch Rejection
DIR_IN  = 'study_directory';
INFOLD  = '4.autorej';
OUTFOLD = '5.visrej'; 
VAR     = 'visrej';

% Inspect participants one at a time
PARTICIPANT = 82; 

% Participants
load(fullfile(DIR_IN,'participants.mat'))

% Load data
file        = participants(PARTICIPANT).name;
filename    = fullfile(DIR_IN,INFOLD,file);
load(filename);

block_index = find(~strcmp({autorej.block},'Practice'));
 
for a = 1:length(block_index)
    
    % Select and remove BAD trials and channels using visual inspection 
    cfg                 = [];
    cfg.method          = 'trial';
    cfg.alim            = 50;
    cfg.eogscale        = 0.1;
    cfg.ecgscale        = 0.05;
    cfg.preproc.detrend = 'yes';
    cfg.preproc.demean  = 'yes';
    visrej(a)           = ft_rejectvisual(cfg,autorej(block_index(a)));
end

% Save output
fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
save(fullpath1,VAR);
    
%% Summary Statistics Trial/Epoch Rejection
DIR_IN  = 'study_directory';
INFOLD  = '5.visrej';
OUTFOLD = '6.allrej'; 
VAR     = 'allrej';

% Inspect participants one at a time
PARTICIPANT = 82; 

% Participants
load(fullfile(DIR_IN,'participants.mat'))

% Load data
file        = participants(PARTICIPANT).name;
filename    = fullfile(DIR_IN,INFOLD,file);
load(filename);

for a = 1:length(visrej)

    % Select and remove BAD trials and channels using summary statistics
    cfg         = [];
    cfg.channel	= {'all','-ECG','-HEOG','-VEOG','-Marker'};
    cfg.method  = 'summary';
    cfg.metric  = 'var'; % or 'zvalue'
    allrej(a)   = ft_rejectvisual(cfg,visrej(a));
end

% Save output
fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
save(fullpath1,VAR);

%% Channel Interpolation
DIR_IN  = 'study_directory';
INFOLD1 = '3.trials';
INFOLD2 = '6.allrej';
OUTFOLD = '7.interp'; 
VAR     = 'interp';

% Load necessary files
elec = ft_read_sens('easycap-M1.txt'); % 3D electrode layout
load('easycapM1.mat') % electrode layout
load('easycapM1_neighb.mat'); % electrode neighbours
load(fullfile(DIR_IN,'participants.mat')) % participants

for p = 1:length(participants)
    
    % Load data
    file        = participants(p).name;
    filename1   = fullfile(DIR_IN,INFOLD1,file);
    filename2   = fullfile(DIR_IN,INFOLD2,file);

    load(filename1);
    load(filename2);

    block_index = find(~strcmp({trials.block},'Practice'));
    trials      = trials(block_index);

    for a = 1:length(allrej)

        % Store list of removed trials and channels
        bad_trial{a} = ismember(trials(a).sampleinfo(:,1), setdiff(trials(a).sampleinfo(:,1), allrej(a).sampleinfo(:,1)));
        bad_chans{a} = setdiff(trials(a).label(1:64),allrej(a).label);

        % Interpolate data for any removed channels
        cfg.method          = 'weighted'; 
        cfg.trials          = 'all';
        cfg.neighbours      = neighbours;
        cfg.elec            = elec;
        cfg.lambda          = 1e-5;
        cfg.order           = 4;
        cfg.missingchannel  = bad_chans{a};

        interp(a) = ft_channelrepair(cfg, allrej(a));

        % Store list of removed trials and channels in final data structure
        interp(a).hdr.rmv_trls = bad_trial{a};
        interp(a).hdr.rmv_chns = bad_chans{a};

        % Add back ECG/EOG/Marker channels from trials, minus bad trials
        interp(a).label(65:68) = trials(a).label(65:68);     
        good_trial = find(~bad_trial{a});

        for trl = 1:length(good_trial)
            interp(a).trial{trl}(65:68,:) = trials(a).trial{good_trial(trl)}(65:68,:);
        end
    end
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% ICA Using Concatenated Structure
DIR_IN  = 'study_directory';
INFOLD  = '7.interp';
OUTFOLD = '8.comps'; 
VAR     = 'comps';

% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    
    % Load data
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    load(filename);

    % Concatenate data structure to allow ICA
    ICA_data    = MDDvsCTRL_ICAdata(interp);

    % ICA
    cfg                 = [];
    cfg.method          = 'runica';
    cfg.channel         = 1:64; % EEG channels only
    
    comps = ft_componentanalysis(cfg,ICA_data);

    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% ICA Component Inspection
DIR_IN  = 'study_directory';
INFOLD  = '8.comps';
OUTFOLD = '8.comps'; 
VAR     = 'comps';

PARTICIPANT = 82; % inspect participants one at a time

% Participants
load(fullfile(DIR_IN,'participants.mat'))

% Load data
file        = participants(PARTICIPANT).name;
filename    = fullfile(DIR_IN,INFOLD,file);
load(filename);

% Inspect ICA components
comps = MDDvsCTRL_ICAcheck(comps);

% Save output
fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
save(fullpath1,VAR);

%% ICA Component Rejection
DIR_IN  = 'study_directory';
INFOLD1 = '7.interp';
INFOLD2 = '8.comps';
OUTFOLD = '9.postica'; 
VAR     = 'postica';

% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    
    % Load data
    file        = participants(p).name;
    filename1   = fullfile(DIR_IN,INFOLD1,file);
    load(filename1);
    
    filename2   = fullfile(DIR_IN,INFOLD2,file);
    load(filename2);
    
    for m = 1:length(interp)
        
        % Remove non-EEG channels
        cfg             = [];
        cfg.channel     = {'all','-ECG','-HEOG','-VEOG','-Marker'};
        interp(m)       = ft_selectdata(cfg,interp(m));
        
        % Create component time series using original data
        cfg             = [];
        cfg.unmixing    = comps.unmixing; % NxN unmixing matrix
        cfg.topolabel   = comps.topolabel; % Nx1 cell-array with the channel labels
        cfg.channel     = 1:64;
        comp_orig       = ft_componentanalysis(cfg,interp(m));

        % Original data reconstructed excluding rejected components
        cfg             = [];
        cfg.component   = comps.rejected;
        cfg.channel     = 1:64;
        postica(m)      = ft_rejectcomponent(cfg,comp_orig,interp(m));
    end

    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% Re-referencing to average
DIR_IN  = 'study_directory';
INFOLD  = '9.postica';
OUTFOLD = '10.reref'; 
VAR     = 'reref';

% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    
    % Load data
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    load(filename);
    
    % Rereference
    clear reref
    cfg             = [];
    cfg.channel     = 'all';
    cfg.refchannel  = 'all';
    cfg.reref       = 'yes';
    cfg.refmethod   = 'avg';
    cfg.trials      = 'all';

    for n = 1:length(postica)
        reref(n) = ft_preprocessing(cfg,postica(n));
    end

    % Add event and block-name structures to file
    for o = 1:length(postica)
        reref(o).event = postica(o).event;
        reref(o).block = postica(o).block;
    end
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% Create IAPS ERPs
DIR_IN  = 'study_directory';
INFOLD  = '10.reref';
OUTFOLD = '11.IAPS'; 
VAR     = 'IAPS';
VALENCE = {'NEGATIVE','NEUTRAL'};
BASE    = [-0.2 0.0];
    
% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    
    % Load data
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    load(filename);
    
    % Find IAPS block
    clear IAPS
    IAPS_index 	= find(strcmp({reref.block},'IAPS'));

    % Average IAPS trials
    for v = 1:length(VALENCE)
        tempdata            = reref(IAPS_index);
        trial_index         = find(strcmp({tempdata.event.type},VALENCE{v}));
        tempdata.trial      = tempdata.trial(trial_index);
        tempdata.time       = tempdata.time(trial_index);
        tempdata.sampleinfo = tempdata.sampleinfo(trial_index,:);

        cfg             = [];
        cfg.channel     = 'all';
        cfg.baseline    = BASE;
        tempdata        = ft_timelockbaseline(cfg,tempdata);

        cfg             = [];
        IAPS(v)         = ft_timelockanalysis(cfg,tempdata);
    end
    
    % Label files
    for v = 1:length(VALENCE)
        IAPS(v).block   = VALENCE{v};
    end
    
    % Save output
    fullpath1 = fullfile(DIR_IN,OUTFOLD,file);
    save(fullpath1,VAR);
end

%% IAPS Statistics and Plotting

% Statistics and potting final results
MDDvsCTRL_IAPS_plots

%% Congratulations!
% Well done on making it all the way to the bottom of the script.
% I hope it wasn't too much of discombobulated mess and that you managed, 
% against all the odds, to make some sense of what I've written.
% Feel free to contact me if there is anything that isn't quite clear.
% Best,
%   Steve
%   stevan.nikolin@unsw.edu.au | stevan.nikolin@gmail.com
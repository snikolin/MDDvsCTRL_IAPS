function data = MDDvsCTRL_preprocess(data)

% Parameters
LINE_NOISE  = 50;

% Hide non-EEG channels from the preprocessing step
EEG_chans   = strcmp(data.hdr.chantype,'EEG');
extra_dat   = data.trial{1}(~EEG_chans,:);
extra_lab   = data.label(~EEG_chans);

% Remove line noise
[b,a] = butter(2, ([LINE_NOISE*0.98 LINE_NOISE*1.02])/(data.fsample * 0.5),'stop');
data.trial{1}(EEG_chans,:) = filtfilt(b,a,data.trial{1}(EEG_chans,:)')';

% Filtering
cfg             = [];
cfg.bpfilter    = 'yes';
cfg.bpfiltord   = 2;
cfg.bpfreq      = [0.5 70];
cfg.bsfilter    = 'yes';
cfg.bsfiltord   = 2;
cfg.bsfreq      = [LINE_NOISE*0.98 LINE_NOISE*1.02]; % Second line-noise removal
cfg.demean      = 'yes';
cfg.detrend     = 'yes';
cfg.channel     = {'all', '-ECG','-HEOG', '-VEOG','-Marker'}';
data            = ft_preprocessing(cfg,data);

% Add extra channels back to data
data.trial{1}(~EEG_chans,:) = extra_dat;
data.label(~EEG_chans)      = extra_lab;

%% Toolbox dependencies
% fieldtrip
% cbrewer
% univarScatter
% padcat
% boundedline

%% Generate data for statistics
DIR_IN      = 'study_directory';
INFOLD      = '11.IAPS';
LPP_WINDOW  = [0.4 1.0]; % Preregistered LPP time interval
CHANNELS    = {'Cz','CPz','Pz','CP1','CP2'}; % Preregistered LPP channels
VALENCE     = {'NEGATIVE','NEUTRAL'};

% Participants
load(fullfile(DIR_IN,'participants.mat'))

for p = 1:length(participants)
    
    % Load data
    file        = participants(p).name;
    filename    = fullfile(DIR_IN,INFOLD,file);
    load(filename);
    
    IAPS = rmfield(IAPS,{'var','dof','block'}); 
    
    % Calculate negative - neutral difference
    cfg             = []; 
    cfg.operation   = 'subtract';
    cfg.parameter   = 'avg';  
    IAPS(3) = ft_math(cfg,IAPS(1),IAPS(2));
    
    % Store all participant data together 
    IAPS_negative{p}    = IAPS(1);
    IAPS_neutral{p}     = IAPS(2);
    IAPS_difference{p}  = IAPS(3);
    
    % Find channels and time interval of interest
    LPP_coi     = ismember(IAPS(1).label,CHANNELS);
    LPP_toi     = IAPS(1).time >= LPP_WINDOW(1) & IAPS(1).time <= LPP_WINDOW(2);
    
    % Create IAPS data structure
    IAPS_data(p).name           = participants(p).name;
    IAPS_data(p).group          = participants(p).group;
    IAPS_data(p).LPP_neg_val    = mean(IAPS(1).avg(LPP_coi,LPP_toi),'all');
    IAPS_data(p).LPP_neut_val   = mean(IAPS(2).avg(LPP_coi,LPP_toi),'all');
    IAPS_data(p).LPP_diff_val   = mean(IAPS(3).avg(LPP_coi,LPP_toi),'all');
    IAPS_data(p).LPP_neg        = mean(IAPS(1).avg(LPP_coi,:),1);
    IAPS_data(p).LPP_neut       = mean(IAPS(2).avg(LPP_coi,:),1);
    IAPS_data(p).LPP_diff       = mean(IAPS(3).avg(LPP_coi,:),1);
end

% Save IAPS_data
fullpath1 = fullfile(DIR_IN,'IAPS_data');
save(fullpath1,'IAPS_data');

%% Calculate Grand Averages
DIR_IN      = 'study_directory';

% Load group list
load(fullfile(DIR_IN,'IAPS_data.mat'))
i_MDD   = [IAPS_data.group] == 0 & [IAPS_data.outlier] ~= 1; % 0 = MDD
i_CTRL  = [IAPS_data.group] == 1 & [IAPS_data.outlier] ~= 1; % 1 = Control

% Set configuration
cfg             = [];
cfg.method      = 'across';
cfg.channel     = 'all';
cfg.latency     = 'all';
cfg.parameter   = 'avg';

% Create grand averages (no outliers)
IAPS_grand_MDD_negative       = ft_timelockgrandaverage(cfg,IAPS_negative{i_MDD});
IAPS_grand_MDD_neutral        = ft_timelockgrandaverage(cfg,IAPS_neutral{i_MDD});
IAPS_grand_MDD_difference     = ft_timelockgrandaverage(cfg,IAPS_difference{i_MDD});
IAPS_grand_control_negative   = ft_timelockgrandaverage(cfg,IAPS_negative{i_CTRL});
IAPS_grand_control_neutral    = ft_timelockgrandaverage(cfg,IAPS_neutral{i_CTRL});
IAPS_grand_control_difference = ft_timelockgrandaverage(cfg,IAPS_difference{i_CTRL});

% Save data
save('IAPS_grandavg','IAPS_grand_MDD_neutral','IAPS_grand_MDD_negative','IAPS_grand_MDD_difference',...
    'IAPS_grand_control_neutral','IAPS_grand_control_negative','IAPS_grand_control_difference')
    
%% Statistics for LPPs
DIR_IN          = 'study_directory';
PERMUTATIONS    = 3000;
ALPHA           = 0.05;
TIMEOFINTEREST  = [0 1];
MIN_NEIGHB      = 2;
TAIL            = 0;

% Load files
load(fullfile(DIR_IN,'participants.mat')); % participants
load('easycapM1.mat') % electrode layout
load('easycapM1_neighb.mat'); % electrode neighbours
load('IAPS_all_data.mat'); % all participant IAPS data
load('elec_easycapM1.mat'); % from: elec = ft_read_sens('easycap-M1.txt')
load('IAPS_data.mat');

% Index of participants
i_MDD   = [IAPS_data.group] == 0 & [IAPS_data.outlier] ~= 1; % 0 = MDD
i_CTRL  = [IAPS_data.group] == 1 & [IAPS_data.outlier] ~= 1; % 1 = Control

% Configuration
cfg                     = [];
cfg.layout              = lay;
cfg.neighbours          = neighbours;
cfg.channel             = 'all'; 
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_indepsamplesT'; 
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum';
cfg.clusterthreshold    = 'nonparametric_individual';
cfg.minnbchan           = MIN_NEIGHB;
cfg.tail                = TAIL;
cfg.clustertail         = TAIL;
cfg.alpha               = ALPHA;
cfg.clusteralpha        = ALPHA;
cfg.numrandomization    = PERMUTATIONS;
cfg.latency             = TIMEOFINTEREST;
cfg.ivar                = 1;
cfg.design              = [ones(1,sum(i_MDD)),2*ones(1,sum(i_CTRL))];

% Negative: MDD vs Healthy
[stats_IAPS_negative]   = ft_timelockstatistics(cfg,IAPS_negative{i_MDD},IAPS_negative{i_CTRL});

% Neutral: MDD vs Healthy
[stats_IAPS_neutral]    = ft_timelockstatistics(cfg,IAPS_neutral{i_MDD},IAPS_neutral{i_CTRL});

% Difference (negative - neutral): MDD vs Healthy
[stats_IAPS_difference] = ft_timelockstatistics(cfg,IAPS_difference{i_MDD},IAPS_difference{i_CTRL});

%% Generating Data - Exploratory Cluster Values
% Identify channels and time intervals from signficant permutation testing and enter below
% These should be obtained from: stats_IAPS_difference
DIR_IN          = 'study_directory';
CLUSTER_WINDOW  = [0.6621 0.7549]; % significant interval from permutation testing
CLUSTER_CHANS   = {'POz','Pz','P1','P2','CPz','CP1'}; % significant channels from permutation testing

for p = 1:length(IAPS_data)
    
    % Identify time and channel indices
    time_index = IAPS_neutral{1}.time >= CLUSTER_WINDOW(1) & IAPS_neutral{1}.time <= CLUSTER_WINDOW(2);
    chan_index = ismember(IAPS_neutral{1}.label,CLUSTER_CHANS);

    % Create IAPS data structure
    IAPS_data(p).LPP_neg_cluster    = mean(IAPS_negative{p}.avg(chan_index,time_index),'all');
    IAPS_data(p).LPP_neut_cluster   = mean(IAPS_neutral{p}.avg(chan_index,time_index),'all');
    IAPS_data(p).LPP_diff_cluster   = mean(IAPS_difference{p}.avg(chan_index,time_index),'all');
end

% Save updated IAPS_data
fullpath1 = fullfile(DIR_IN,'IAPS_data');
save(fullpath1,'IAPS_data');

%% ERPs - Contrast Unpleasant vs Neutral Valences - Preregistered Cluster
% Figures 2A, 2B, and 2C

load('IAPS_time.mat')

CMAP        = cbrewer('qual','Set1', 3, 'cubic'); % Unpleasant | [] | Neutral
HEIGHT      = 4;
WIDTH       = 8;
FONTSIZE    = 8;
YADJ        = 1.01;
TIME_WINDOW = [-0.2 1.0];
BOOTSTRAPS  = 1000;
LPP_POS     = 0.8;
ERP_TIME    = find(IAPS_time >= TIME_WINDOW(1) & IAPS_time <= TIME_WINDOW(2));

% Generate data
for n = 1:length(IAPS_data)
    LPP_trace_neg(n,:)    = IAPS_data(n).LPP_neg(ERP_TIME);
    LPP_trace_neut(n,:)   = IAPS_data(n).LPP_neut(ERP_TIME);
    LPP_trace_diff(n,:)   = IAPS_data(n).LPP_diff(ERP_TIME);
end

% Index of participants
i_MDD   = [IAPS_data.group] == 0 & [IAPS_data.outlier] ~= 1; % 0 = MDD
i_CTRL  = [IAPS_data.group] == 1 & [IAPS_data.outlier] ~= 1; % 1 = Control

% MDD
temp1 = []; temp2 = []; error1 = []; error2 = [];
temp1 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_neg(i_MDD,:));
temp2 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_neut(i_MDD,:));
temp1 = sort(temp1);
temp2 = sort(temp2);
error1(:,1)  = temp1(975,:) - temp1(500,:);
error1(:,2)  = temp1(500,:) - temp1(25,:); 
error2(:,1)  = temp2(975,:) - temp2(500,:);
error2(:,2)  = temp2(500,:) - temp2(25,:);

figure
h1 = boundedline(IAPS_time(ERP_TIME),temp1(500,:),error1,...
    'cmap',CMAP(1,:),'alpha','transparency',0.2,'linewidth',1.5);
h2 = boundedline(IAPS_time(ERP_TIME),temp2(500,:),error2,...
    'cmap',CMAP(3,:),'alpha','transparency',0.2,'linewidth',1.5);
set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
[~,L] = legend({'Unpleasant','Neutral'}, 'location','southeast');
PatchInLegend = findobj(L,'type','patch');
set(PatchInLegend(1), 'FaceAlpha', 1);
set(PatchInLegend(2), 'FaceAlpha', 1);

hold on
yerp_max = max([max(temp1),max(temp2)])*YADJ;
yerp_min = min([min(temp1),min(temp2)])*YADJ;
plot([IAPS_time(ERP_TIME(1)) IAPS_time(ERP_TIME(end))],...
    [0 0],'LineStyle','--','Color','black','HandleVisibility','off');
plot([0 0],[yerp_min yerp_max],'LineStyle','--','Color','black','HandleVisibility','off');
rectangle('Position', [0.4 yerp_min 0.6 (yerp_max-yerp_min)],...
    'FaceColor', [0.1 0.1 0.1 0.1],'EdgeColor','none') 

axis tight
xlabel('Time (secs)')
ylabel('Amplitude (uV)')
title('MDD')
text(0.7,yerp_max*LPP_POS,'LPP','FontWeight','bold','Fontsize',FONTSIZE,'VerticalAlignment','bottom')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT])
print -dtiff LPP_trace_MDDvalence_hypcluster.tif -r300

% Controls
temp1 = []; temp2 = []; error1 = []; error2 = [];
temp1 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_neg(i_CTRL,:));
temp2 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_neut(i_CTRL,:));
temp1 = sort(temp1);
temp2 = sort(temp2);
error1(:,1)  = temp1(975,:) - temp1(500,:);
error1(:,2)  = temp1(500,:) - temp1(25,:); 
error2(:,1)  = temp2(975,:) - temp2(500,:);
error2(:,2)  = temp2(500,:) - temp2(25,:);

figure
h1 = boundedline(IAPS_time(ERP_TIME),temp1(500,:),error1,...
    'cmap',CMAP(1,:),'alpha','transparency',0.2,'linewidth',1.5);
h2 = boundedline(IAPS_time(ERP_TIME),temp2(500,:),error2,...
    'cmap',CMAP(3,:),'alpha','transparency',0.2,'linewidth',1.5);
set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
[~,L] = legend({'Unpleasant','Neutral'}, 'location','southeast');
PatchInLegend = findobj(L,'type','patch');
set(PatchInLegend(1), 'FaceAlpha', 1);
set(PatchInLegend(2), 'FaceAlpha', 1);

hold on
yerp_max = max([max(temp1),max(temp2)])*YADJ;
yerp_min = min([min(temp1),min(temp2)])*YADJ;
plot([IAPS_time(ERP_TIME(1)) IAPS_time(ERP_TIME(end))],...
    [0 0],'LineStyle','--','Color','black','HandleVisibility','off');
plot([0 0],[yerp_min yerp_max],'LineStyle','--','Color','black','HandleVisibility','off');
rectangle('Position', [0.4 yerp_min 0.6 (yerp_max-yerp_min)],...
    'FaceColor', [0.1 0.1 0.1 0.1],'EdgeColor','none') 

axis tight
xlabel('Time (secs)')
ylabel('Amplitude (uV)')
title('Controls')
text(0.7,yerp_max*LPP_POS,'LPP','FontWeight','bold','Fontsize',FONTSIZE,'VerticalAlignment','bottom')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT])
print -dtiff LPP_trace_CTRLvalence_hypcluster.tif -r300

% Difference
CMAP        = cbrewer('qual','Dark2', 3, 'cubic'); % MDD | [] | Control

temp1 = []; temp2 = []; error1 = []; error2 = [];
temp1 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_diff(i_MDD,:));
temp2 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_diff(i_CTRL,:));
temp1 = sort(temp1);
temp2 = sort(temp2);
error1(:,1)  = temp1(975,:) - temp1(500,:);
error1(:,2)  = temp1(500,:) - temp1(25,:); 
error2(:,1)  = temp2(975,:) - temp2(500,:);
error2(:,2)  = temp2(500,:) - temp2(25,:);

figure
h1 = boundedline(IAPS_time(ERP_TIME),temp1(500,:),error1,...
    'cmap',CMAP(1,:),'alpha','transparency',0.2,'linewidth',1.5);
h2 = boundedline(IAPS_time(ERP_TIME),temp2(500,:),error2,...
    'cmap',CMAP(3,:),'alpha','transparency',0.2,'linewidth',1.5);
set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
[~,L] = legend({'MDD','Control'}, 'location','southeast');
PatchInLegend = findobj(L,'type','patch');
set(PatchInLegend(1), 'FaceAlpha', 1);
set(PatchInLegend(2), 'FaceAlpha', 1);

hold on
yerp_max = max([max(temp1),max(temp2)])*YADJ;
yerp_min = min([min(temp1),min(temp2)])*YADJ;
plot([IAPS_time(ERP_TIME(1)) IAPS_time(ERP_TIME(end))],...
    [0 0],'LineStyle','--','Color','black','HandleVisibility','off');
plot([0 0],[yerp_min yerp_max],'LineStyle','--','Color','black','HandleVisibility','off');
rectangle('Position', [0.4 yerp_min 0.6 (yerp_max-yerp_min)],...
    'FaceColor', [0.1 0.1 0.1 0.1],'EdgeColor','none') 

axis tight
xlabel('Time (secs)')
ylabel('Amplitude (uV)')
title('Difference: unpleasant - neutral')
text(0.7,yerp_max*LPP_POS,'LPP','FontWeight','bold','Fontsize',FONTSIZE,'VerticalAlignment','bottom')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT])
print -dtiff LPP_trace_difference_hypcluster.tif -r300

%% ERPs - Contrast Unpleasant vs Neutral Valences - Exploratory Cluster
% Figures 3A, 3B, and 3C

load('IAPS_time.mat')

CMAP        = cbrewer('qual','Set1', 3, 'cubic'); % Unpleasant | [] | Neutral
HEIGHT      = 4;
WIDTH       = 8;
FONTSIZE    = 8;
YADJ        = 1.01;
TIME_WINDOW = [-0.2 1.0];
BOOTSTRAPS  = 1000;
LPP_POS     = 0.8;
ERP_TIME    = find(IAPS_time >= TIME_WINDOW(1) & IAPS_time <= TIME_WINDOW(2));
CLUSTER_CH  = {'POz','Pz','P1','P2','CPz','CP1'};
CLUSTER_T   = [0.6621 0.7549];

% Find time and channel
time_index  = IAPS_neutral{1}.time >= CLUSTER_T(1) & IAPS_neutral{1}.time <= CLUSTER_T(2);
chan_index  = ismember(IAPS_neutral{1}.label,CLUSTER_CH);

% Generate data
for n = 1:length(IAPS_data)
    LPP_trace_neg(n,:)    = mean(IAPS_negative{n}.avg(chan_index,ERP_TIME),1);
    LPP_trace_neut(n,:)   = mean(IAPS_neutral{n}.avg(chan_index,ERP_TIME),1);
    LPP_trace_diff(n,:)   = mean(IAPS_difference{n}.avg(chan_index,ERP_TIME),1);
end

% Index of participants
i_MDD   = [IAPS_data.group] == 0 & [IAPS_data.outlier] ~= 1; % 0 = MDD
i_CTRL  = [IAPS_data.group] == 1 & [IAPS_data.outlier] ~= 1; % 1 = Control

% MDD
temp1 = []; temp2 = []; error1 = []; error2 = [];
temp1 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_neg(i_MDD,:));
temp2 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_neut(i_MDD,:));
temp1 = sort(temp1);
temp2 = sort(temp2);
error1(:,1)  = temp1(975,:) - temp1(500,:);
error1(:,2)  = temp1(500,:) - temp1(25,:); 
error2(:,1)  = temp2(975,:) - temp2(500,:);
error2(:,2)  = temp2(500,:) - temp2(25,:);

figure
h1 = boundedline(IAPS_time(ERP_TIME),temp1(500,:),error1,...
    'cmap',CMAP(1,:),'alpha','transparency',0.2,'linewidth',1.5);
h2 = boundedline(IAPS_time(ERP_TIME),temp2(500,:),error2,...
    'cmap',CMAP(3,:),'alpha','transparency',0.2,'linewidth',1.5);
set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
[~,L] = legend({'Unpleasant','Neutral'}, 'location','northeast');
PatchInLegend = findobj(L,'type','patch');
set(PatchInLegend(1), 'FaceAlpha', 1);
set(PatchInLegend(2), 'FaceAlpha', 1);

hold on
yerp_max = max([max(temp1),max(temp2)])*YADJ;
yerp_min = min([min(temp1),min(temp2)])*YADJ;
plot([IAPS_time(ERP_TIME(1)) IAPS_time(ERP_TIME(end))],...
    [0 0],'LineStyle','--','Color','black','HandleVisibility','off');
plot([0 0],[yerp_min yerp_max],'LineStyle','--','Color','black','HandleVisibility','off');
rectangle('Position', [CLUSTER_T(1) yerp_min (CLUSTER_T(2) - CLUSTER_T(1)) (yerp_max-yerp_min)],...
    'FaceColor', [0.1 0.1 0.1 0.1],'EdgeColor','none') 

axis tight
xlabel('Time (secs)')
ylabel('Amplitude (uV)')
title('MDD')
text(mean(CLUSTER_T),yerp_min+0.5,'LPP','FontWeight','bold',...
    'Fontsize',FONTSIZE,'VerticalAlignment','bottom','HorizontalAlignment','center')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT])
print -dtiff LPP_trace_MDDvalence_sigcluster.tif -r300

% Controls
temp1 = []; temp2 = []; error1 = []; error2 = [];
temp1 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_neg(i_CTRL,:));
temp2 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_neut(i_CTRL,:));
temp1 = sort(temp1);
temp2 = sort(temp2);
error1(:,1)  = temp1(975,:) - temp1(500,:);
error1(:,2)  = temp1(500,:) - temp1(25,:); 
error2(:,1)  = temp2(975,:) - temp2(500,:);
error2(:,2)  = temp2(500,:) - temp2(25,:);

figure
h1 = boundedline(IAPS_time(ERP_TIME),temp1(500,:),error1,...
    'cmap',CMAP(1,:),'alpha','transparency',0.2,'linewidth',1.5);
h2 = boundedline(IAPS_time(ERP_TIME),temp2(500,:),error2,...
    'cmap',CMAP(3,:),'alpha','transparency',0.2,'linewidth',1.5);
set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
[~,L] = legend({'Unpleasant','Neutral'}, 'location','northeast');
PatchInLegend = findobj(L,'type','patch');
set(PatchInLegend(1), 'FaceAlpha', 1);
set(PatchInLegend(2), 'FaceAlpha', 1);

hold on
yerp_max = max([max(temp1),max(temp2)])*YADJ;
yerp_min = min([min(temp1),min(temp2)])*YADJ;
plot([IAPS_time(ERP_TIME(1)) IAPS_time(ERP_TIME(end))],...
    [0 0],'LineStyle','--','Color','black','HandleVisibility','off');
plot([0 0],[yerp_min yerp_max],'LineStyle','--','Color','black','HandleVisibility','off');
rectangle('Position', [CLUSTER_T(1) yerp_min (CLUSTER_T(2) - CLUSTER_T(1)) (yerp_max-yerp_min)],...
    'FaceColor', [0.1 0.1 0.1 0.1],'EdgeColor','none') 

axis tight
xlabel('Time (secs)')
ylabel('Amplitude (uV)')
title('Controls')
text(mean(CLUSTER_T),yerp_min+0.5,'LPP','FontWeight','bold',...
    'Fontsize',FONTSIZE,'VerticalAlignment','bottom','HorizontalAlignment','center')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT])
print -dtiff LPP_trace_CTRLvalence_sigcluster.tif -r300

% Difference
CMAP        = cbrewer('qual','Dark2', 3, 'cubic'); % MDD | [] | Control

temp1 = []; temp2 = []; error1 = []; error2 = [];
temp1 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_diff(i_MDD,:));
temp2 = bootstrp(BOOTSTRAPS, @(x) mean(x), LPP_trace_diff(i_CTRL,:));
temp1 = sort(temp1);
temp2 = sort(temp2);
error1(:,1)  = temp1(975,:) - temp1(500,:);
error1(:,2)  = temp1(500,:) - temp1(25,:); 
error2(:,1)  = temp2(975,:) - temp2(500,:);
error2(:,2)  = temp2(500,:) - temp2(25,:);

figure
h1 = boundedline(IAPS_time(ERP_TIME),temp1(500,:),error1,...
    'cmap',CMAP(1,:),'alpha','transparency',0.2,'linewidth',1.5);
h2 = boundedline(IAPS_time(ERP_TIME),temp2(500,:),error2,...
    'cmap',CMAP(3,:),'alpha','transparency',0.2,'linewidth',1.5);
set(get(get(h1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
[~,L] = legend({'MDD','Control'}, 'location','northwest');
PatchInLegend = findobj(L,'type','patch');
set(PatchInLegend(1), 'FaceAlpha', 1);
set(PatchInLegend(2), 'FaceAlpha', 1);

hold on
yerp_max = max([max(temp1),max(temp2)])*YADJ;
yerp_min = min([min(temp1),min(temp2)])*YADJ;
plot([IAPS_time(ERP_TIME(1)) IAPS_time(ERP_TIME(end))],...
    [0 0],'LineStyle','--','Color','black','HandleVisibility','off');
plot([0 0],[yerp_min yerp_max],'LineStyle','--','Color','black','HandleVisibility','off');
rectangle('Position', [CLUSTER_T(1) yerp_min (CLUSTER_T(2) - CLUSTER_T(1)) (yerp_max-yerp_min)],...
    'FaceColor', [0.1 0.1 0.1 0.1],'EdgeColor','none') 

axis tight
xlabel('Time (secs)')
ylabel('Amplitude (uV)')
title('Difference: unpleasant - neutral')
text(mean(CLUSTER_T),yerp_max*LPP_POS,'LPP','FontWeight','bold',...
    'Fontsize',FONTSIZE,'VerticalAlignment','bottom','HorizontalAlignment','center')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT])
print -dtiff LPP_trace_difference_sigcluster.tif -r300

%% TOPOGRAPHY - LPP for Preregistered Cluster
% Figures 2A, 2B, 2C, and 2E

LPP_WINDOW  = [0.4 1];
CHANNELS    = {'Cz','CPz','Pz','CP1','CP2'};
HEIGHT      = 4;

% Generate group difference values
cfg                         = [];
cfg.operation               = 'subtract';
cfg.parameter               = 'avg';
IAPS_grand_DIFF_difference  = ft_math(cfg,IAPS_grand_MDD_difference,...
    IAPS_grand_control_difference);

% Configuration details
cfg                     = [];
cfg.xlim                = LPP_WINDOW;
cfg.layout              = 'easycapM1.mat';
cfg.comment             = 'xlim';
cfg.commentpos          = 'title';
cfg.highlight           = 'on';
cfg.highlightchannel    = CHANNELS;
cfg.highlightcolor      = 'red';
cfg.highlightsymbol     = '*';
cfg.highlightsize       = 6;
cfg.highlightfontsize   = 6;

% Group Valence Values
cfg.zlim                = [-2.2 2.25];

figure; ft_topoplotER(cfg,IAPS_grand_MDD_negative); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_MDD_neg_hypcluster.tif -r300

figure; ft_topoplotER(cfg,IAPS_grand_MDD_neutral); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_MDD_neut_hypcluster.tif -r300

figure; ft_topoplotER(cfg,IAPS_grand_control_negative); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_CTRL_neg_hypcluster.tif -r300

figure; ft_topoplotER(cfg,IAPS_grand_control_neutral); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_CTRL_neut_hypcluster.tif -r300

% Valence Differences
cfg.zlim                = [-0.75 0.5];

figure; ft_topoplotER(cfg,IAPS_grand_MDD_difference);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_MDD_diff_hypcluster.tif -r300

figure; ft_topoplotER(cfg,IAPS_grand_control_difference); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_CTRL_diff_hypcluster.tif -r300

% Valence and Group Difference
cfg.zlim                = [-0.45 0.6];

ft_topoplotER(cfg,IAPS_grand_DIFF_difference); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_DIFF_diff_hypcluster.tif -r300

%% TOPOGRAPHY - LPP for Exploratory Cluster
% Figures 3A, 3B, 3C, and 3E

% Find significant channels 
for c = 1:length(stats_IAPS_difference.label)
    cluster_channels(c,:) = stats_IAPS_difference.negclusterslabelmat(c,:) == 1;
end

cluster_channels        = sum(cluster_channels,2);
[A,B]                   = sort(cluster_channels,'descend');
cluster_sig_channels    = stats_IAPS_difference.label(B);

% Find significant time
for t = 1:length(stats_IAPS_difference.time)
    cluster_time(t,:) = stats_IAPS_difference.negclusterslabelmat(:,t) == 1;
end
cluster_time            = sum(cluster_time,2);
cluster_time_index      = find(logical(cluster_time));
cluster_time_coords     = [cluster_time_index(1) cluster_time_index(end)];

% Parameters
HEIGHT          = 4;
CLUSTER_WINDOW  = stats_IAPS_difference.time(cluster_time_coords); % 0.6621 0.7549

% Configuration details
cfg                     = [];
cfg.xlim                = CLUSTER_WINDOW;
cfg.layout              = 'easycapM1.mat';
cfg.comment             = 'xlim';
cfg.commentpos          = 'title';
cfg.highlight           = 'on';
cfg.highlightchannel    = cluster_sig_channels(1:6);
cfg.highlightcolor      = 'red';
cfg.highlightsymbol     = '*';
cfg.highlightsize       = 6;
cfg.highlightfontsize   = 6;

% Group Valence Values
cfg.zlim                = [-2.05 1.85];

figure; ft_topoplotER(cfg,IAPS_grand_MDD_negative); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_MDD_neg_sigcluster.tif -r300

figure; ft_topoplotER(cfg,IAPS_grand_MDD_neutral); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_MDD_neut_sigcluster.tif -r300

figure; ft_topoplotER(cfg,IAPS_grand_control_negative); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_CTRL_neg_sigcluster.tif -r300

figure; ft_topoplotER(cfg,IAPS_grand_control_neutral);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_CTRL_neut_sigcluster.tif -r300

% Valence Differences
cfg.zlim                = [-1.15 0.9];

figure; ft_topoplotER(cfg,IAPS_grand_MDD_difference);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_MDD_diff_sigcluster.tif -r300

figure; ft_topoplotER(cfg,IAPS_grand_control_difference); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_CTRL_diff_sigcluster.tif -r300

% Valence and Group Difference
cfg.zlim                = [-2.5 2];

ft_topoplotER(cfg,IAPS_grand_DIFF_difference); 
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 1.2*HEIGHT 1.2*HEIGHT])
print -dtiff LPP_topo_DIFF_diff_sigcluster.tif -r300

%% SCATTERPLOTs - Contrast Unpleasant vs Neutral Valences - Preregistered Cluster
% Figures 2D

CMAP        = cbrewer('qual','Set1', 3, 'cubic'); % Unpleasant | [] | Neutral
colours     = vertcat(CMAP(1,:),CMAP(3,:));
COMPRESS    = 7;
HEIGHT      = 4;
WIDTH       = 6;
FONTSIZE    = 8;
SCATTERSZ   = 5;
MARKWIDTH   = 0.6;
YADJ        = 1.05;
XLIMS       = [0.5 2.5];
XLIMDIFF    = [0.5 1.5];
LINEW       = 0.5;
YLIMITS     = [-1.6547 1.9108];

% Load IAPS data
load('IAPS_data.mat')

% Generate figures for preregistered cluster
i           = [IAPS_data.group] == 0; % 0 is MDD; 1 = Control
LPP_neg     = [IAPS_data.LPP_neg_val]; LPP_neg([IAPS_data.outlier]) = NaN;
LPP_neut    = [IAPS_data.LPP_neut_val]; LPP_neut([IAPS_data.outlier]) = NaN;
ylimits     = [min([LPP_neg LPP_neut])*YADJ max([LPP_neg LPP_neut])*YADJ];

% MDD plot
figure
UnivarScatter([LPP_neg(i); LPP_neut(i)]','Label',{'U','N'},...
    'MarkerFaceColor',colours,'MarkerEdgeColor','black','Width',MARKWIDTH,...
    'Whiskers','box','Compression',COMPRESS,'PointSize',SCATTERSZ); 
hold on
ylabel('LPP Amplitude (uV)')
xlim(XLIMS)
ylim(ylimits)
title('MDD')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT],'papersize',[21.0 29.7])
print -dtiff IAPS_MDD_hypcluster.tif -r300

% Control plot
figure
UnivarScatter([LPP_neg(~i); LPP_neut(~i)]','Label',{'U','N'},...
    'MarkerFaceColor',colours,'MarkerEdgeColor','black','Width',MARKWIDTH,...
    'Whiskers','box','Compression',COMPRESS,'PointSize',SCATTERSZ); 
hold on
ylabel('LPP Amplitude (uV)')
xlim(XLIMS)
ylim(ylimits)
title('Controls')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT],'papersize',[21.0 29.7])
print -dtiff IAPS_Control_hypcluster.tif -r300s

%% SCATTERPLOTs - Contrast Unpleasant vs Neutral Valences - Exploratory Cluster
% Figures 3D

CMAP        = cbrewer('qual','Set1', 3, 'cubic'); % Unpleasant | [] | Neutral
colours     = vertcat(CMAP(1,:),CMAP(3,:));
COMPRESS    = 7;
HEIGHT      = 4;
WIDTH       = 6;
FONTSIZE    = 8;
SCATTERSZ   = 5;
MARKWIDTH   = 0.6;
YADJ        = 1.05;
XLIMS       = [0.5 2.5];
XLIMDIFF    = [0.5 1.5];
LINEW       = 0.5;
YLIMITS     = [-1.6547 1.9108];

% Load IAPS data
load('IAPS_data.mat')

% Generate figures for exploratory cluster
i           = [IAPS_data.group] == 0; % 0 is MDD; 1 = Control
LPP_neg     = [IAPS_data.LPP_neg_cluster]; LPP_neg([IAPS_data.outlier]) = NaN;
LPP_neut    = [IAPS_data.LPP_neut_cluster]; LPP_neut([IAPS_data.outlier]) = NaN;
ylimits     = [min([LPP_neg LPP_neut])*YADJ max([LPP_neg LPP_neut])*YADJ];

% MDD plot
figure
UnivarScatter([LPP_neg(i); LPP_neut(i)]','Label',{'U','N'},...
    'MarkerFaceColor',colours,'MarkerEdgeColor','black','Width',MARKWIDTH,...
    'Whiskers','box','Compression',COMPRESS,'PointSize',SCATTERSZ); 
hold on
ylabel('LPP Amplitude (uV)')
xlim(XLIMS)
ylim(ylimits)
title('MDD')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT],'papersize',[21.0 29.7])
print -dtiff IAPS_MDD_sigcluster.tif -r300

% Control plot
figure
UnivarScatter([LPP_neg(~i); LPP_neut(~i)]','Label',{'U','N'},...
    'MarkerFaceColor',colours,'MarkerEdgeColor','black','Width',MARKWIDTH,...
    'Whiskers','box','Compression',COMPRESS,'PointSize',SCATTERSZ); 
hold on
ylabel('LPP Amplitude (uV)')
xlim(XLIMS)
ylim(ylimits)
title('Controls')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT],'papersize',[21.0 29.7])
print -dtiff IAPS_Control_sigcluster.tif -r300s

%% SCATTERPLOTs - Emotional Regulation Strategy Plot
CMAP        = cbrewer('qual','Dark2', 3, 'cubic'); % MDD | [] | Control
colours     = vertcat(CMAP(1,:),CMAP(3,:));
COMPRESS    = 7;
HEIGHT      = 4;
WIDTH       = 6;
FONTSIZE    = 8;
SCATTERSZ   = 5;
MARKWIDTH   = 0.6;
YADJ        = 1.05;
XLIMS       = [0.5 2.5];
XLIMDIFF    = [0.5 1.5];
LINEW       = 0.5;
YLIMITS     = [-1.6547 1.9108];

% Load IAPS data
load('IAPS_data.mat')

% Generate figure
LPP_neg         = [IAPS_data.LPP_neg_val]; LPP_neg([IAPS_data.outlier]) = NaN;

% No strategy index
i_MDD           = [IAPS_data.group] == 0 & [IAPS_data.strategy] == 0; % 0 is MDD; 1 = Control
i_CTRL          = [IAPS_data.group] == 1 & [IAPS_data.strategy] == 0; % 0 is No Strategy; 1 = Strategy

% No strategy plot
ylimits = [min(LPP_neg)*YADJ max(LPP_neg)*YADJ];

figure
UnivarScatter(padcat(LPP_neg(i_MDD), LPP_neg(i_CTRL))','Label',{'MDD','Control'},...
    'MarkerFaceColor',colours,'MarkerEdgeColor','black','Width',MARKWIDTH,...
    'Whiskers','box','Compression',COMPRESS,'PointSize',SCATTERSZ); 
hold on
ylabel('LPP Amplitude (uV)')
xlim(XLIMS)
ylim(ylimits)
title('No strategy')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT],'papersize',[21.0 29.7])
print -dtiff IAPS_nostrategy.tif -r300

% Strategy index
i_MDD           = [IAPS_data.group] == 0 & [IAPS_data.strategy] == 1; % 0 is MDD; 1 = Control
i_CTRL          = [IAPS_data.group] == 1 & [IAPS_data.strategy] == 1; % 0 is No Strategy; 1 = Strategy

% Strategy plot
ylimits = [min(LPP_neg)*YADJ max(LPP_neg)*YADJ];

figure
UnivarScatter(padcat(LPP_neg(i_MDD), LPP_neg(i_CTRL))','Label',{'MDD','Control'},...
    'MarkerFaceColor',colours,'MarkerEdgeColor','black','Width',MARKWIDTH,...
    'Whiskers','box','Compression',COMPRESS,'PointSize',SCATTERSZ); 
hold on
ylabel('LPP Amplitude (uV)')
xlim(XLIMS)
ylim(ylimits)
title('Strategy')
set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT],'papersize',[21.0 29.7])
print -dtiff IAPS_strategy.tif -r300

%% SCATTERPLOTs - Anxiety Plot
CMAP        = cbrewer('qual','Dark2', 3, 'cubic'); % MDD | [] | Control
HEIGHT      = 4;
WIDTH       = 4;
FONTSIZE    = 8;
SCATTER_SZ  = 5; % Size of scatter dots
SCATTER_THK = 0.5; % Thickness of scatter dots
SHAPE       = 'o'; % Shape of vertex points (e.g. 'x')

% Index
i           = [IAPS_data.group] == 0; % 0 is MDD; 1 = Control

LPP_neg         = [IAPS_data.LPP_neg_val]; LPP_neg([IAPS_data.outlier]) = NaN;
anxiety         = [IAPS_data.DASS_anx];

% Anxiety plot
scatter(anxiety(i), LPP_neg(i),...
        SCATTER_SZ,CMAP(1,:),'filled',SHAPE,'Linewidth',SCATTER_THK,...
    'MarkerEdgeColor','black');

temp = [anxiety(i); LPP_neg(i)]';
temp(any(isnan(temp),2),:) = [];

P  = polyfit(temp(:,1), temp(:,2),1);
x0 = min(temp(:,1)); 
x1 = max(temp(:,1));
xi = linspace(x0,x1) ;
yi = P(1)*xi+P(2);
hold on
plot(xi,yi,'linestyle',':','color','black','linewidth',1) ;

box on

ylabel('LPP Amplitude (uV)')
xlabel('DASS Anxiety')
xlim([0 32])

set(gca,'FontSize',FONTSIZE)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 WIDTH HEIGHT],'papersize',[21.0 29.7])
print -dtiff IAPS_anxiety.tif -r300

function comps = MDDvsCTRL_ICAcheck(comps)
% Identify which components are to be removed
% Done through automated process to suggest components for removal 
% with follow-up visual confirmation
% Plot, inspect, and identify first 25 ICs for removal (Display 5 x 5)
% At the prompt, list the components to reject
% Example: [3,6,7,8,9];

%% Automated ICA inspection
[blinks, lateral, muscle] = MDDvsCTRL_autoICA(comps);
suggested_comps = [blinks, lateral, muscle];

%% Visual ICA inspection
MDDvsCTRL_ICAplot(comps,suggested_comps);

%% List components for rejection
% Suggestion of components to remove
disp('Automated ICA reject components are:')
disp(['Blinks: ', num2str(blinks)]);
disp(['Lateral eye movements: ', num2str(lateral)]);
disp(['Muscle: ', num2str(muscle)]);

% Creates prompt for user to include which components need to be rejected
prompt          = 'list components to reject. Approximately 10-15% of total:';
x               = input(prompt);
comps.rejected  = x;
close

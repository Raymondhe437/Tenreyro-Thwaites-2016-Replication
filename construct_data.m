%POC to make sure we can modify and track changes to MATLAB via VS Code and then upload the modifications onto Github
%POC number 2 to make it work without needing VS Code

startS=89;  %
endS=224;
samplelength=endS-startS+1; % Sample period - dates of shocks - need to ensure data runs L periods earlier and H periods later

% Import the data
[~, ~, raw] = xlsread('data.xlsx','Policy','A2:A255');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
FFR = reshape([raw{:}],size(raw));

%% Clear temporary variables
clearvars raw R;

[~, ~, raw] = xlsread('data.xlsx','Shocks','A2:B255');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
E = reshape([raw{:}],size(raw));

%% Clear temporary variables
clearvars raw R;

[~, ~, raw] = xlsread('data.xlsx','alt_states','A2:D255');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
alt_states = reshape([raw{:}],size(raw));

%% Clear temporary variables
clearvars raw R;

%% Import the data
[~, ~, raw] = xlsread('data.xlsx','LHS','A2:AM255');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
LHS = reshape([raw{:}],size(raw));
dates = LHS(:,1);
LHS = LHS(:,2:end);

%% Clear temporary variables
clearvars raw R;

R=FFR;                      % Policy variable being instrumented
[~, ~, names] = xlsread('data.xlsx','LHS','B1:AM1');
names(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),names)) = {''};

save data.mat %R E LHS alt_states names dates startS endS samplelength


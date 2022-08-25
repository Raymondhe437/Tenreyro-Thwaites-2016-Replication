%% This is the file with annotation. -hsk-
%% Github
startS=89;  %
endS=224;
samplelength=endS-startS+1; % Sample period - dates of shocks - need to ensure data runs L periods earlier and H periods later

% Import the data
[~, ~, raw] = xlsread('C:\Users\HUSANG KIM\Desktop\Papers\3. MPU&MPP\6. Local_Projection_MATLAB\data.xlsx','Policy','A2:A255');
%hsk% ~ seems like a filler to fill arguments inside [ ]
%hsk% [num, txt, raw] = xlsread(filename, sheet, cellrange)
%hsk% num: only stores numeric values from the excel file
%hsk% txt: only stores text values from the excel file
%hsk% raw: stores every value as it is

raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
%hsk% The policy sheet is empty : may be this constructs an empty
%hsk% cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw) shows the
%result of applying the function @(x) to each cell of raw. 
%hsk& & and && are both AND operator. && employs short-circuiting behaviour.
%hsk% find cells that are not empty, numeric, and nan simultaneously. I
%guess nans are also treated as numeric.

%dataset(R or Python's NULL equivalent)
%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
FFR = reshape([raw{:}],size(raw));
%hsk% raw returns the 254 x 1 array
%hsk% raw{:}: prints out each value in each array as if running a for loop
%hsk% [raw{:}]: lists each array in Column through Column.
%hsk% size(raw) returns the size of array which is 254 x 1
%hsk% reshape([raw{:}],size(raw)) creates a vector of output variable.
%% Clear temporary variables
clearvars raw R;

[~, ~, raw] = xlsread('C:\Users\HUSANG KIM\Desktop\Papers\3. MPU&MPP\6. Local_Projection_MATLAB\data.xlsx','Shocks','A2:B255');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
E = reshape([raw{:}],size(raw));

%% Clear temporary variables
clearvars raw R;

[~, ~, raw] = xlsread('C:\Users\HUSANG KIM\Desktop\Papers\3. MPU&MPP\6. Local_Projection_MATLAB\data.xlsx','alt_states','A2:D255');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
alt_states = reshape([raw{:}],size(raw));

%% Clear temporary variables
clearvars raw R;

%% Import the data
[~, ~, raw] = xlsread('C:\Users\HUSANG KIM\Desktop\Papers\3. MPU&MPP\6. Local_Projection_MATLAB\data.xlsx','LHS','A2:AM255');
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
[~, ~, names] = xlsread('C:\Users\HUSANG KIM\Desktop\Papers\3. MPU&MPP\6. Local_Projection_MATLAB\data.xlsx','LHS','B1:AM1');
names(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),names)) = {''};

save data.mat %R E LHS alt_states names dates startS endS samplelength

%hsk%
% R = FFR: Federal funds rates (policy rates)
% E : Shocks
% LHS: macro variables
% alt_states: 4 different alternative state variables.
% names: names of LHS macro variables
% dates: dates
% StartS endS samplelength


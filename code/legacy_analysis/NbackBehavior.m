% HISTORICAL REFERENCE ONLY. Requires original inputs/toolboxes and path configuration.
% Use code/RUN_ALL.m for the portable release analyses.
% Load the Excel file
filename = 'E:\Mito_DICOM\Event\MiSBIE Eprime Old Computer 10-26-23.xlsx';
% filename = 'E:\Mito_DICOM\Event\MiSBIE Eprime New Computer 10-26-23.xlsx';
data = readtable(filename);

% Extract relevant columns
participantID = data{:, 'Subject'};  % Assuming 'B' is the second column
targetResponses = data{:, 'Slide13_ACC'};  % Assuming 'KO' is the 291st column

% Identify relevant rows (assuming n-back task related rows are marked)
% You may need to update this condition based on how n-back rows are identified
nbackRows = ~isnan(targetResponses);

% Filter out non-n-back task rows
participantID = participantID(nbackRows);
targetResponses = targetResponses(nbackRows);

% Determine the number of participants
uniqueParticipants = unique(participantID);
numParticipants = length(uniqueParticipants);

% Initialize an array to store the accuracy for each participant
accuracy = zeros(numParticipants, 1);

% Compute accuracy for each participant
for i = 1:numParticipants
    % Get current participant ID
    currentParticipantID = uniqueParticipants(i);
    
    % Extract trials for the current participant
    participantTrials = targetResponses(participantID == currentParticipantID);
    
    % Compute the number of blocks and trials
    numTrials = length(participantTrials);
    numBlocks = numTrials / 10;  % Assuming 10 trials per block
    
    % Calculate the accuracy
    correctResponses = sum(participantTrials);
    accuracy(i) = correctResponses / numTrials;
end

% Create a table to store the results
results = table(uniqueParticipants, accuracy, 'VariableNames', {'ParticipantID', 'Accuracy'});

% Display the results
disp(results);


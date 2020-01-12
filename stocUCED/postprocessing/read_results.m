function [res, output] = read_results(contentDir)

files = dir(fullfile(char(contentDir)));

% clean up weird files
idx = [];
for f=1:length(files)
    if (strcmp(files(f).name,'.') || ... 
        strcmp(files(f).name,'..') || ... 
        strcmp(files(f).name,'.DS_Store') || ...
        strcmp(files(f).name,'notes.txt') || ... 
        strcmp(files(f).name(end-3:end),'.zip'))
        idx = [idx; f];
    end
end
files(idx) = [];

renIdx = (236:266); % renewable generators
hydIdx = (173:215);
conIdx = (1:235);   % conventional generators
STGenIdx = [101:171, 216:222];  % ST-UC generators
DAGenIdx = [1:100, 172:215, 223:266]; % DA-UC generators

% read rest of them
for f=1:length(files)
   res.names{f} = files(f).name;   
   res.stats{f} = read_stats( strcat(contentDir, files(f).name, '/') );      
   res.gen{f} = read_result( strcat(contentDir, files(f).name, '/'), '_genED' );
   % renCurtailment contains ONLY SOLAR + WIND GENS
   res.renCurtailment{f} = read_result_by_idx(strcat(contentDir, files(f).name, '/'), '_overGenED', renIdx);
   % conOverGen ONLY CONVENTIONAL GENS
   res.conOverGen{f} = read_result_by_idx( strcat(contentDir, files(f).name, '/'), '_overGenED', conIdx);
   % contains all generators
   res.extraGen{f} = read_result(strcat(contentDir, files(f).name, '/'), '_overGenED');
   res.allCommitments{f} = read_result( strcat(contentDir, files(f).name, '/'), '_commitments');
   res.DAGenCommitments{f} = read_result_by_idx( strcat(contentDir, files(f).name, '/'), '_commitments', DAGenIdx);
   res.STGenCommitments{f} = read_result_by_idx( strcat(contentDir, files(f).name, '/'), '_commitments', STGenIdx);
   res.renNhydroGen{f} = read_result_by_idx(strcat(contentDir, files(f).name, '/'), '_genED', [hydIdx, renIdx]);
   res.renNhydroGen{f} = res.renNhydroGen{f} - [res.extraGen{f}(hydIdx,:); res.renCurtailment{f}];
end

% postprocess
[pwCost, pwCO2, pwNOX, pwSO2] = getPiecewiseLinearCostAndGHGEstimates(res);
for f=1:length(files)   
   res.stats{f}.convOverGen = sum(res.conOverGen{f})';
   res.stats{f}.renCurtailment = sum(res.renCurtailment{f})';
   res.stats{f}.percentSTGenCommit = (sum(res.STGenCommitments{f})./(sum(res.DAGenCommitments{f})+sum(res.STGenCommitments{f})))';
   res.stats{f}.percentSTGen = (sum(res.gen{f}(STGenIdx,:)) ./ sum(res.gen{f}))';
   res.stats{f}.renNhydroGen = sum(res.renNhydroGen{f})';
   res.stats{f}.renNhydroPercent = (res.stats{f}.renNhydroGen - res.stats{f}.renCurtailment) ./ ((sum(res.gen{f}) - sum(res.extraGen{f}))');   
   
   res.stats{f}.pwTotalCost = sum(pwCost{f})' + res.stats{f}.NoLoadCost + res.stats{f}.StartUpCost;
   res.stats{f}.pwCO2 = sum(pwCO2{f})';
   res.stats{f}.pwNOX = sum(pwNOX{f})';
   res.stats{f}.pwSO2 = sum(pwSO2{f})';   
end

% % % write the stats
for f=1:length(files)
    writetable(struct2table(res.stats{f}), 'subhourly_summarized_outputs.xlsx', 'FileType', 'spreadsheet', 'Sheet', res.names{f})
end

end

function res = read_result(contentDir, type)
% Reads all solution files,
% puts them in a cell array

% read files
% contentDir = "/results";
files = dir(fullfile(char(strcat(contentDir, '*', type, '.sol'))));

% allocate memory
res = [];

% read file contents
for f=1:length(files)
    res = [res, csvread(strcat(contentDir,files(f).name))];
end
end

function res = read_result_by_idx(contentDir, type, idx)
% Reads all solution files,
% puts them in a cell array

% read files
% contentDir = "/results";
files = dir(fullfile(char(strcat(contentDir, '*', type, '.sol'))));

% allocate memory
res = [];

% read file contents
for f=1:length(files)
    tmp = csvread(strcat(contentDir,files(f).name));
    tmp = tmp(idx,:);
    res = [res, tmp];
end
end

function res = read_stats(contentDir)
% Reads all stat solution files,
% puts them in a cell array

% read files
% contentDir = "/results";
files = dir(fullfile(char(strcat(contentDir, '*stats.sol'))));

% read file contents
res = tdfread(strcat(contentDir,files(1).name));
for f=2:length(files)
    tdf = tdfread(strcat(contentDir,files(f).name));
    fields = fieldnames(tdf);
    for i = 1:numel(fields)
        res.(fields{i}) = [res.(fields{i}); tdf.(fields{i})];
    end
end

v = repmat(1:length(files),[96 1]);
res.day = v(:);

end
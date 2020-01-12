function printGenerationMixes(res, casename)
% Example Usage:
% printGenerationMixes(res, 'Ren3.00_Res0.10_0.025')

% read generator fuel data
fid = fopen('NREL118_GeneratorFuels.csv');
data = textscan(fid, '%s%s', 'delimiter', ',', 'headerlines', 1);
fclose(fid);

for f=1:length(res.names)
    if strcmp(res.names{f}(5:end), casename)
        printGenerationMix(res, f, data);
    end
end
end

function printGenerationMix(res, idx, data)

% initialize inputs and basic parameters
numGen = length(data{1});
genIdx = zeros(numGen, 1);
for i=1:numGen
   genIdx(i) = str2double(data{1}{i});
end

raw_gen = res.gen{idx};
raw_overgen = res.extraGen{idx};

fuelTypes = unique(data{2});
numFuels = length(fuelTypes);

numPeriods = size(raw_gen);
numPeriods = numPeriods(2);

Gen = zeros(numFuels, numPeriods);
OverGen = zeros(numFuels, numPeriods);

for p=1:numPeriods
    for f=1:numFuels
        Gen(f, p) = sum(raw_gen(strcmp(data{2},fuelTypes(f)), p));
        OverGen(f, p) = sum(raw_overgen(strcmp(data{2},fuelTypes(f)), p));
    end
end

output = [];
numDays = 7;
v = repmat(1:numDays,[numPeriods/numDays 1]);
output.days = v(:);
output.numPeriods = (1:numPeriods)';
for f=1:numFuels
    fuelTypes{f} = regexprep(fuelTypes{f}, ' ', '_');
    output = setfield(output, fuelTypes{f}, (Gen(f,:) - OverGen(f,:))');
end
output.convOverGen = zeros(numPeriods, 1);
output.swCurtail = zeros(numPeriods, 1);

for f=1:numFuels
    if (strcmp(fuelTypes{f}, 'Solar') || strcmp(fuelTypes{f}, 'Wind'))
        output.swCurtail = output.swCurtail + OverGen(f,:)';
    else 
        output.convOverGen = output.convOverGen + OverGen(f,:)';
    end
end

% create the table and print it
writetable(struct2table(output), 'genMix.xlsx', 'FileType', 'spreadsheet', 'Sheet', res.names{idx})
end
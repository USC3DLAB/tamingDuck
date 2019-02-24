function [pwCost, pwCO2, pwNOX, pwSO2] = getPiecewiseLinearCostAndGHGEstimates(res)

% get basic parameters
nFiles = length(res.stats);
numGen = size(res.gen{1});
numGen = numGen(1);
numBPs = 5;

% get piecewise linear function coefficients 
[PCC, CO2, NOX, SO2, BPs] = loadCoefs(numGen, numBPs);

pwCost = cell(nFiles, 1);
pwCO2 = cell(nFiles, 1);
pwNOX = cell(nFiles, 1);
pwSO2 = cell(nFiles, 1);

for j=1:nFiles
    [pwCost{j}, pwCO2{j}, pwNOX{j}, pwSO2{j}] = getPiecewiseLinearCostAndGHGEstimate(res.gen{j}, res.allCommitments{j}, 0.25, PCC, CO2, NOX, SO2, BPs);
end

display('Warning: Make sure cost-statistics include no-load (intercept) and start up costs. They must be postprocessed');

end

function [pwCost, pwCO2, pwNOX, pwSO2] = getPiecewiseLinearCostAndGHGEstimate (generation, commitment, hoursPerPeriod, PCC, CO2, NOX, SO2, BPs)

% basic params
dim = size(generation);
numGen = dim(1);
numPeriods = dim(2);
numBPs = 5;
EPS = 1e-6;

% compute the outputs %
pwCost = zeros(numGen, numPeriods);
pwCO2 = zeros(numGen, numPeriods);
pwNOX = zeros(numGen, numPeriods);
pwSO2 = zeros(numGen, numPeriods);

% intercept (no-load GHG emissions; no-load costs to be added later)
pwCO2 = pwCO2 + commitment .* CO2(:,1);
pwNOX = pwNOX + commitment .* NOX(:,1);
pwSO2 = pwSO2 + commitment .* SO2(:,1);

% remaining linear pieces
for g=1:numGen
    for t=1:numPeriods
        remGen = generation(g,t);
        for b=1:numBPs
            if remGen > EPS 
               if b==1
                   val = min(remGen, BPs(g,b));
               else
                   val = min(remGen, BPs(g,b)-BPs(g,b-1));
               end
               
               remGen = remGen - val;

               pwCost(g,t) = pwCost(g,t) + val * PCC(g,b+1);
               pwCO2(g,t) = pwCO2(g,t) + val * CO2(g,b+1);
               pwNOX(g,t) = pwNOX(g,t) + val * NOX(g,b+1);
               pwSO2(g,t) = pwSO2(g,t) + val * SO2(g,b+1);
            else
                break
            end
        end
    end
end 

pwCost = pwCost * hoursPerPeriod;
pwCO2 = pwCO2 * hoursPerPeriod;
pwNOX = pwNOX * hoursPerPeriod;
pwSO2 = pwSO2 * hoursPerPeriod;

end

function [PCC, CO2, NOX, SO2, BPs] = loadCoefs(numGen, numBPs)

% parse data
heatdata = importdata('NREL118_GeneratorHeatRates.csv', ',');
fueldata = importdata('NREL118_Fuels and emission rates.csv', ',');

% create a map to read emission rate indices
dim         = size(fueldata.textdata);
keySet      = cell(dim(1)-1,1);
valueSet    = zeros(dim(1)-1,1);
for c=2:dim(1)
    %{fueldata.textdata{c,1}, c}
    keySet{c-1} = fueldata.textdata{c,1};
    valueSet(c-1) = c-1;
end
fuel2idx = containers.Map(keySet, valueSet);

% read data
PCC = zeros(numGen, numBPs+1); % piecewise cost coefficients
CO2 = zeros(numGen, numBPs+1);
NOX = zeros(numGen, numBPs+1);
SO2 = zeros(numGen, numBPs+1);
BPs = zeros(numGen, numBPs); % BreakPointS of the piecewise linear functions

convfactor = 0.453592;

for g=1:numGen
    if (~isKey(fuel2idx, heatdata.textdata{g+1,2})) 
        continue;
    end
    
    b=1;
    PCC(g,b) = heatdata.data(g,b)*fueldata.data(fuel2idx(heatdata.textdata{g+1,2}),1);
    CO2(g,b) = heatdata.data(g,b)*fueldata.data(fuel2idx(heatdata.textdata{g+1,2}),2)*convfactor;
    NOX(g,b) = heatdata.data(g,b)*fueldata.data(fuel2idx(heatdata.textdata{g+1,2}),3)*convfactor;
    SO2(g,b) = heatdata.data(g,b)*fueldata.data(fuel2idx(heatdata.textdata{g+1,2}),4)*convfactor;
    
    for b=2:(numBPs+1)
        PCC(g,b) = heatdata.data(g,b)*fueldata.data(fuel2idx(heatdata.textdata{g+1,2}),1)/1000;
        CO2(g,b) = heatdata.data(g,b)*fueldata.data(fuel2idx(heatdata.textdata{g+1,2}),2)/1000*convfactor;
        NOX(g,b) = heatdata.data(g,b)*fueldata.data(fuel2idx(heatdata.textdata{g+1,2}),3)/1000*convfactor;
        SO2(g,b) = heatdata.data(g,b)*fueldata.data(fuel2idx(heatdata.textdata{g+1,2}),4)/1000*convfactor;
    end
end

dim = size(heatdata.data); 
BPs(1:dim(1),:) = heatdata.data(:,(numBPs+2):(numBPs*2+1)); % these two may not be of the same size due to blank GHG/cost coefficients on the bottom of the CSV

% additional V&OM charge ($/MWh)
PCC(1:dim(1),2:end) = PCC(1:dim(1),2:end) + heatdata.data(:,end);

% remove NAs
CO2(isnan(CO2)) = 0;
NOX(isnan(NOX)) = 0;
SO2(isnan(SO2)) = 0;
PCC(isnan(PCC)) = 0;
BPs(isnan(BPs)) = 99999;
end
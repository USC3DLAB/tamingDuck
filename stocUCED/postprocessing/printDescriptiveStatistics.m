function output = printDescriptiveStatistics(res)    
    % cost 
    writeDescriptiveStatistics(res, 'pwTotalCost', 'pwTotalCost');
    
    % commitment-related
    writeDescriptiveStatistics(res, 'x0x23ofOnGen', 'NumOnlineGen');
    writeDescriptiveStatistics(res, 'percentSTGenCommit', '%OnlineSTGen');
    
    % generation-related
    writeDescriptiveStatistics(res, 'UsedGen', 'UsedGen');
    writeDescriptiveStatistics(res, 'LoadShed', 'UnmetLoad');
    writeDescriptiveStatistics(res, 'renCurtailment', 'SWCurtail');
    writeDescriptiveStatistics(res, 'convOverGen', 'OverGen');
    writeDescriptiveStatistics(res, 'renNhydroGen', 'SWHGen');
    writeDescriptiveStatistics(res, 'renNhydroPercent', '%SWH');
    
    % environmental impacts
    writeDescriptiveStatistics(res, 'pwCO2', 'CO2');
    writeDescriptiveStatistics(res, 'pwNOX', 'NOX');
    writeDescriptiveStatistics(res, 'pwSO2', 'SO2');
    
    % system-characteristics
    writeDescriptiveStatistics(res, 'MinGenReq', 'MinGenReq');
    writeDescriptiveStatistics(res, 'RampUpCap', 'RampUpCap');
    writeDescriptiveStatistics(res, 'RampDownCap', 'RampDownCap');
    
end

function output = writeDescriptiveStatistics(res, field, oname)
    output.names = res.names';    
    output.ops = [];
    output.SW_coef = [];
    output.UC_res_coef = [];
    output.ED_res_coef = [];
    output.Battery_coef = [];
    
    output.mean = [];
    output.std = [];
    output.min = [];
    output.max = [];
    for f=1:length(res.stats)
        instance = strsplit(output.names{f}, '_');
        output.ops = [output.ops; instance{1}];
        output.SW_coef = [output.SW_coef; str2double(instance{2}(4:end))];
        output.UC_res_coef = [output.UC_res_coef; str2double(instance{3}(4:end))];
        output.ED_res_coef = [output.ED_res_coef; str2double(instance{4})];
        if isempty(strfind(output.names{f},'Bt'))
            output.Battery_coef = [output.Battery_coef; 0.0];
        else
            output.Battery_coef = [output.Battery_coef; str2num(output.names{f}(end))];
        end
        
        output.min = [output.min; min(getfield(res.stats{f}, field))];
        output.max = [output.max; max(getfield(res.stats{f}, field))];
        output.mean = [output.mean; mean(getfield(res.stats{f}, field))];
        output.std = [output.std; std(getfield(res.stats{f}, field))];
    end
    
    writetable(struct2table(output), 'statistics_15minResolution.xlsx', 'FileType', 'spreadsheet', 'Sheet', oname)
end
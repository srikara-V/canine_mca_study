file = "C:\RSI 2023\project work\data sorted or external\Filtered_canine_illumina_calls_pheno (1).csv";

data = readtable(file, 'Format', 'auto');

types = unique(data.chrom);

for j=1:size(types, 1)

    rowsToKeep = strcmp(data.chrom, types(j));
    dataOfSample = data(rowsToKeep, :);
    scatterData = zeros(size(dataOfSample, 1), 4);
    
    for i=1:size(dataOfSample,1)
        health = 0;
        if contains(dataOfSample.HEALTH_STATUS(i), 'LSA') || contains(dataOfSample.HEALTH_STATUS(i), 'OSA') || contains(dataOfSample.HEALTH_STATUS(i), 'sarcoma') || contains(dataOfSample.HEALTH_STATUS(i), 'ymphoma') || contains(dataOfSample.HEALTH_STATUS(i), 'MCT') 
            health = 1;
        elseif strcmp(dataOfSample.HEALTH_STATUS(i), 'NA')
            health = 0.5;
        else
            health = 0;
        end
        
        scatterData(i,:) = [dataOfSample.beg_canFam3_txt(i), dataOfSample.end_canFam3_txt(i), dataOfSample.length(i), health];
    end

    figure;
    scatter(scatterData(:,1), scatterData(:,2), [], scatterData(:,4), 'filled');
    % Add a colorbar for the labels
    cbar = colorbar;
    cbar.Label.String = 'Label';
    % Set labels for each axis
    xlabel('Start of mCA');
    ylabel('End of mCA');
    title(types(j));

    saveas(gcf, "scatter_plot"+types(j)+".png");
end


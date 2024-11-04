file = 'C:\RSI 2023\project work\Dog_MoChA\canine_illumina.calls.tsv';
data = readtable(file, 'FileType', 'text', 'Delimiter', '\t');

samples = unique(data.sample_id);
chromosomes = unique(data.chrom);
types = unique(data.type);

counts = [];
totalCounts = zeros(size(chromosomes, 1), size(types, 1));
for i=1:size(chromosomes, 1)
    for j=1:size(types, 1)
        rowsToKeepX = strcmp(data.type, types{j});
        x = data(rowsToKeepX, :);

        rowsToKeepY = strcmp(x.chrom, chromosomes{i});
        y = x(rowsToKeepY, :);
        counts = horzcat(counts, size(y,1));

        x = [];
        
    end
    totalCounts(i,:) = counts;
    counts = [];
end

for i=1:size(types,1)
    figure;
    bar(totalCounts(:,i));
    title(types{i});
end
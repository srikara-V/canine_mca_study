%% TODO: Replace path
file1 = "C:\Davidson\canine_illumina.sample_phenotypes.txt";
file2 = "C:\Davidson\canine_illumina.calls.filtered_original_IDs.tsv";

data1 = readtable(file1, "FileType", "text");
data2 = readtable(file2, 'FileType', 'text');

Samples = data1.mocha_sample_id;
Cancer = data1.Cancer;

chromosomes = unique(data2.chrom);
types = unique(data2.type);

for i=1:size(chromosomes,1)
    filter = data2(strcmp(data2.chrom, chromosomes(i)),:);
    for j=1:size(types,1)
        filter = filter(strcmp(filter.type, types(j)),:);

        has_mCAs = [];
        for k=1:size(Samples,1)
            if size(filter(strcmp(filter.sample_id, Samples(k)),:),1) ~= 0
                has_mCAs = [has_mCAs; {'Y'}];
            else
                has_mCAs = [has_mCAs; {'N'}];
            end
        end

        t = table(Samples, Cancer, has_mCAs);
        writetable(t, "table_"+chromosomes(i)+"_"+types(j)+".csv");

        filter = data2(strcmp(data2.chrom, chromosomes(i)),:);
    end
end
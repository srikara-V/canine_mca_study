%% TODO: Replace path
file1 = "C:\Davidson\canine_illumina.sample_phenotypes.txt";
file2 = "C:\Davidson\canine_illumina.calls.filtered_original_IDs.tsv";

data1 = readtable(file1, "FileType", "text");
data2 = readtable(file2, 'FileType', 'text');

types = {'CN-LOH', 'Gain', 'Loss', 'Undetermined'};

permutations = 5000;

p_values = 0;

%% TODO: Replace path
file = "C:\Davidson\table_chr8_CN-LOH.csv";
data = readtable(file, "Format", "auto");

in_box = data(strcmp(data.has_mCAs, "Y"),:);

if size(in_box) ~= 0
    cancer_in_box = in_box(strcmp(in_box.Cancer, "Y"),:);
    
    p1 = size(cancer_in_box,1) / size(in_box,1);
    E_cancer_given_mca = -p1*log2(p1)-(1-p1)*log2(1-p1);

    if p1 == 1 && size(in_box,1) > 4
        E_cancer_given_mca = 0;
    elseif p1==1
        E_cancer_given_mca = 1;
    end
    
    if p1 == 0 && size(in_box,1) > 4
        E_cancer_given_mca = 0;
    elseif p1 == 0
        E_cancer_given_mca = 1;
    end
    
    fake_entropies = zeros([permutations,1]);
    fake_probs = zeros([permutations,1]);

    for k=1:permutations
        fake_box = data(randi([1,4134],1,size(in_box,1)),:);
        cancer_in_fake_box = fake_box(strcmp(fake_box.Cancer, "Y"),:);

        p2 = size(cancer_in_fake_box,1) / size(in_box,1);
        E_cancer_given_mca_fake = -p2*log2(p2)-(1-p2)*log2(1-p2);

        if p2 == 1 && size(in_box,1) > 4
            E_cancer_given_mca_fake = 0;
        elseif p1==1
            E_cancer_given_mca_fake = 1;
        end
        
        if p2 == 0 && size(in_box,1) > 4
            E_cancer_given_mca_fake = 0;
        elseif p2 == 0
            E_cancer_given_mca_fake = 1;
        end

        fake_entropies(k) = E_cancer_given_mca_fake;
        fake_probs(k) = p2;
    end
    p_values = ranksum(E_cancer_given_mca, fake_entropies);
end
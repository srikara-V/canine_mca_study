entropies = zeros([38,4]);
entropies_in = zeros([38,4]);
entropies_out = zeros([38,4]);

types = {'CN-LOH', 'Gain', 'Loss', 'Undetermined'};

for i = 1:38
    for j=1:4
        %% TODO: Replace path
        file = "C:\Davidson\table_chr"+i+"_"+types{j}+".csv";
        data = readtable(file, "Format", "auto");

        in_box = data(strcmp(data.has_mCAs, "Y"),:);
        cancer_in_box = in_box(strcmp(in_box.Cancer, "Y"),:);
        out_box = setxor(data, in_box);
        cancer_out_box = out_box(strcmp(out_box.Cancer, "Y"),:);

        p = size(data(strcmp(data.Cancer, "Y"),:),1) / size(data,1);
        E_cancer = -p*log2(p)-(1-p)*log2(1-p);

        if p == 1 && size(cancer_in_box,1) > 1
            E_cancer = 0;
        end
        
        p1 = size(cancer_in_box,1) / size(in_box,1);
        E_cancer_given_mca = -p1*log2(p1)-(1-p1)*log2(1-p1);

        p2 = size(cancer_out_box,1) / size(out_box,1);
        E_cancer_given_not_mca = -p2*log2(p2)-(1-p2)*log2(1-p2);

        if p1 == 1 && size(cancer_in_box,1) > 4
            E_cancer_given_mca = 0;
        elseif p1==1
            E_cancer_given_mca = E_cancer;
        end

        if p2 == 1 && size(cancer_in_box,1) > 4
            E_cancer_given_not_mca = 0;
        elseif p2==1
            E_cancer_given_not_mca = E_cancer;
        end

        E = E_cancer - E_cancer_given_mca;

        entropies(i,j) = E;
        entropies_in(i,j) = E_cancer_given_mca;
        entropies_out(i,j) = E_cancer_given_not_mca;
    end
end

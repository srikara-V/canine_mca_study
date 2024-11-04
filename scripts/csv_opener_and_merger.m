files = readlines('GEO\GEO_IDs.txt');

c = readtable('GEO\GSE49357_df_1.csv', 'Delimiter', ',');

x = {};
for i=2:size(files,1)
    x{i} = readtable("GEO\"+files(i), 'Delimiter', ',').characteristics_ch1;
end
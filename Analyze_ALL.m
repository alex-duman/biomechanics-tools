% Run All Traces
Plot = 'N'; % set to 'Y' or 'Yes' to have all figures output, otherwise set to any other string or integer

% Read in a list of all trials
list = readtable('List of Reflex Trials.xlsx');

rows = 6:size(list,1);
for i = rows
    T = Analysis_and_Plot_func(char(list{i,1}),Plot);
    if i == min(rows)
        Table = T;
    else
        Table = [Table; T];
    end
end

% writetable(Table, 'Analyzed Reflex Trials.csv')
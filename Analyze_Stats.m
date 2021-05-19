%% Analyze_Stats
% This file will do a repeated measures ANOVA to analyze the stretch reflex
% data recorded throughout the experiment.

%% Reading in Data
D = readtable('Analyzed Reflex Trials.csv');
D.Condition = categorical(D.Condition);
D.Individual = categorical(D.Individual);

%% Performing Summary Calculations
Conditions = categories(D.Condition);
Conditions = [Conditions(end-1:end);Conditions(1:4)]; % rearrange conditions to chronological order
Individuals = categories(D.Individual);

% Determining average values for each individual at each condition (remove
% pseudoreplications) prior to calculating addional statistics
for i = 1:length(Individuals)
    D_Indv = D(D.Individual == Individuals{i},:);
    for j = 1:length(Conditions)
        D_subsamp = D_Indv(D_Indv.Condition==Conditions{j},:);
        if size(D_subsamp,1)>=1
            Individual = D_subsamp.Individual(1);
            Condition = D_subsamp.Condition(1);
            Mean_I_flex = mean(D_subsamp.I_flexion,'omitnan');
            Mean_I_ext = mean(D_subsamp.I_extension,'omitnan');
            Ratio_I_e_f = mean(D_subsamp.Ratio_I_ext_flex(isfinite(D_subsamp.Ratio_I_ext_flex)),'omitnan');
            t_on_flex = mean(D_subsamp.t_on_flex,'omitnan');
            t_on_ext = mean(D_subsamp.t_on_ext,'omitnan');
            Ratio_t_on_e_f = mean(D_subsamp.Ratio_t_on_ext_flex(isfinite(D_subsamp.Ratio_t_on_ext_flex)),'omitnan');
            
            % Creating Summary Table for Individuals
            if i==1 && j==1
                Indv_Means = table(Individual,Condition,Mean_I_flex,Mean_I_ext,Ratio_I_e_f,t_on_flex,t_on_ext,Ratio_t_on_e_f);
            else
                Indv_Means = [Indv_Means; table(Individual,Condition,Mean_I_flex,Mean_I_ext,Ratio_I_e_f,t_on_flex,t_on_ext,Ratio_t_on_e_f)];
            end
        end
    end
end

% Determining Condition means (after removing pseudoreplicants - extra
% trials for each individual at the same condition).
for i = 1:length(Conditions)
    D_subsamp = Indv_Means(Indv_Means.Condition==Conditions{i},:);
    
    % Condition
    Condition = D_subsamp.Condition(1);
    
    % Calculating Averages
    Mean_I_flex = mean(D_subsamp.Mean_I_flex,'omitnan');
    Mean_I_ext = mean(D_subsamp.Mean_I_ext,'omitnan');
    Ratio_I_e_f = mean(D_subsamp.Ratio_I_e_f(isfinite(D_subsamp.Ratio_I_e_f)),'omitnan');
    t_on_flex = mean(D_subsamp.t_on_flex,'omitnan');
    t_on_ext = mean(D_subsamp.t_on_ext,'omitnan');
    Ratio_t_on_e_f = mean(D_subsamp.Ratio_t_on_e_f(isfinite(D_subsamp.Ratio_t_on_e_f)),'omitnan');
    
    % Standard Deviation
    SD_I_flex = std(D_subsamp.Mean_I_flex,'omitnan');
    SD_I_ext = std(D_subsamp.Mean_I_ext,'omitnan');
    SD_Ratio_I_e_f = std(D_subsamp.Ratio_I_e_f(isfinite(D_subsamp.Ratio_I_e_f)),'omitnan');
    SD_t_on_flex = std(D_subsamp.t_on_flex,'omitnan');
    SD_t_on_ext = std(D_subsamp.t_on_ext,'omitnan');
    SD_Ratio_t_on_e_f = std(D_subsamp.Ratio_t_on_e_f(isfinite(D_subsamp.Ratio_t_on_e_f)),'omitnan');
    
    % Making Tables
    if i ==1
        Cond_Means = table(Condition,Mean_I_flex,Mean_I_ext,Ratio_I_e_f,t_on_flex,t_on_ext,Ratio_t_on_e_f);
        Cond_SD = table(Condition,SD_I_flex,SD_I_ext,SD_Ratio_I_e_f,SD_t_on_flex,SD_t_on_ext,SD_Ratio_t_on_e_f);
    else
        Cond_Means = [Cond_Means; table(Condition,Mean_I_flex,Mean_I_ext,Ratio_I_e_f,t_on_flex,t_on_ext,Ratio_t_on_e_f)];
        Cond_SD = [Cond_SD; table(Condition,SD_I_flex,SD_I_ext,SD_Ratio_I_e_f,SD_t_on_flex,SD_t_on_ext,SD_Ratio_t_on_e_f)];
    end
end

%% Summary Plots
% Reordering Conditions
Cond_Means.Condition = reordercats(Cond_Means.Condition,{Conditions{1},Conditions{2},Conditions{3},Conditions{4},Conditions{5},Conditions{6}});
Cond_SD.Condition = reordercats(Cond_SD.Condition,{Conditions{1},Conditions{2},Conditions{3},Conditions{4},Conditions{5},Conditions{6}});
D.Condition = reordercats(D.Condition,{Conditions{1},Conditions{2},Conditions{3},Conditions{4},Conditions{5},Conditions{6}});

% Intensity Plots
figure('Name','Intensity Summary Plots');
subplot(3,1,1);
bar(Cond_Means.Condition,Cond_Means.Mean_I_flex);
hold on
plot(D.Condition,D.I_flexion,'r.','MarkerSize',4);
er1 = errorbar(Cond_Means.Condition,Cond_Means.Mean_I_flex,Cond_SD.SD_I_flex,'Color',[0,0,0],'LineStyle','none','CapSize',18,'LineWidth',0.75);
hold off
ylabel({'Flex. Intensity'; '(V*s)'});

subplot(3,1,2);
bar(Cond_Means.Condition,Cond_Means.Mean_I_ext);
hold on
plot(D.Condition,D.I_extension,'r.',"MarkerSize",4);
er2 = errorbar(Cond_Means.Condition,Cond_Means.Mean_I_ext,Cond_SD.SD_I_ext,'Color',[0,0,0],'LineStyle','none','CapSize',18);
hold off
ylabel({'Ext. Intensity';'(V*s)'});

subplot(3,1,3);
bar(Cond_Means.Condition,Cond_Means.Ratio_I_e_f);
hold on
plot(D.Condition,D.Ratio_I_ext_flex,'r.',"MarkerSize",4);
er3 = errorbar(Cond_Means.Condition,Cond_Means.Ratio_I_e_f,Cond_SD.SD_Ratio_I_e_f,'Color',[0,0,0],'LineStyle','none','CapSize',18);
hold off
ylabel({'Intensity Ratio';'(ext/flex)'});


% Activation Time Plots
figure('Name','Activation Timing Summary Plots');
subplot(3,1,1);
bar(Cond_Means.Condition,Cond_Means.t_on_flex);
hold on
plot(D.Condition,D.t_on_flex,'r.','MarkerSize',4);
er1 = errorbar(Cond_Means.Condition,Cond_Means.t_on_flex,Cond_SD.SD_t_on_flex,'Color',[0,0,0],'LineStyle','none','CapSize',18,'LineWidth',0.75);
hold off
ylabel({'Fraction Time Active';'during Flexion';'(s/s)'});

subplot(3,1,2);
bar(Cond_Means.Condition,Cond_Means.t_on_ext);
hold on
plot(D.Condition,D.t_on_ext,'r.',"MarkerSize",4);
er2 = errorbar(Cond_Means.Condition,Cond_Means.t_on_ext,Cond_SD.SD_t_on_ext,'Color',[0,0,0],'LineStyle','none','CapSize',18);
hold off
ylabel({'Fraction Active time';'during Extension';'(s/s)'});

subplot(3,1,3);
bar(Cond_Means.Condition,Cond_Means.Ratio_t_on_e_f);
hold on
plot(D.Condition,D.Ratio_t_on_ext_flex,'r.',"MarkerSize",4);
er3 = errorbar(Cond_Means.Condition,Cond_Means.Ratio_t_on_e_f,Cond_SD.SD_Ratio_t_on_e_f,'Color',[0,0,0],'LineStyle','none','CapSize',18);
hold off
ylabel({'Ratio Active Time';'(ext/flex)'});
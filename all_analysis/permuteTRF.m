% Permutation testing for TRFs.
load("flux_and_envelope_TRF.mat");

%{
Anova: 
row factor is Tempo (fast/slow)
column factor is Performance Medium (Live/Rec)
A = 
    [
    [LiveFast_part1 RecordedFast_part1],
    [LiveFast_part2 RecordedFast_part2],
    ...
    ...
    ...
    [LiveFast_part21 RecordedFast_part21],
    [LiveSlow_part1 RecordedSlow_part1]
    ...
    ...
    ...
    [LiveSlow_part21 RecordedSlow_part21]
    ]
%}
electrodes = 1:29;
p_Interaction_envelope = zeros(length(time_lags),1);
p_Interaction_flux = zeros(length(time_lags),1);

p_Performance_envelope = zeros(length(time_lags),1);
p_Performance_flux = zeros(length(time_lags),1);

p_Tempo_envelope = zeros(length(time_lags),1);
p_Tempo_flux = zeros(length(time_lags),1);
%%
for t = 1:length(time_lags)
    Anova_envelope = zeros(42,2);
    Anova_flux = zeros(42,2);
    envelope_data = squeeze(mean(all_envelope_TRFs(:,:,electrodes,t),3));
    flux_data = squeeze(mean(all_flux_TRFs(:,:,electrodes,t),3)); % >> both should result in a 21x4 matrix
    %LiveFast condition
    Anova_envelope(1:21,1) = envelope_data(:,1);
    Anova_flux(1:21,1) = flux_data(:,1);
    %LiveSlow condition
    Anova_envelope(22:42,1) = envelope_data(:,2);
    Anova_flux(22:42,1) = flux_data(:,2);
    %RecordedFast condition
    Anova_envelope(1:21,2) = envelope_data(:,3);
    Anova_flux(1:21,2) = flux_data(:,3);
    %RecordedSlow condition
    Anova_envelope(22:42,2) = envelope_data(:,4);
    Anova_flux(22:42,2) = flux_data(:,4);

    [~,~,stats_envelope] = anova2(Anova_envelope,21,'off');
    [~,~,stats_flux] = anova2(Anova_flux,21,'off');

    performanceCompare_envelope = multcompare(stats_envelope,"Estimate","column");
    performanceCompare_flux = multcompare(stats_flux,"Estimate","column");
    tempoCompare_envelope = multcompare(stats_envelope,"Estimate","row");
    tempoCompare_flux = multcompare(stats_flux,"Estimate","row");

    p_Interaction_envelope(t) = stats_envelope.pval;
    p_Interaction_flux(t) = stats_flux.pval;
    p_Performance_envelope(t) = performanceCompare_envelope(end);
    p_Performance_flux(t) = performanceCompare_flux(end);
    p_Tempo_envelope(t) = tempoCompare_envelope(end);
    p_Tempo_flux(t) = tempoCompare_flux(end);
end
%%
 region = time_lags(find(p_Tempo_flux < .01)); % Gives times ~ 46.875 - 62.5 ms, ind 7-9


%% Plot


figure;
hold on
plot(time_lags, p_Interaction_flux, 'DisplayName','Interaction Pval','LineWidth',2.5,'Color',[.5 .8 .5])
plot(time_lags, p_Performance_flux, 'DisplayName','Performance Condition Pval','LineWidth',2.5,'Color',[.8 0 .5])
plot(time_lags, p_Tempo_flux, 'DisplayName','Tempo Condition Pval','LineWidth',2.5,'Color',[1 .6 0])

xregion(min(region),max(region),'DisplayName','p<0.01 for Tempo')

title("Flux Pvals across Time Lags") 
xlabel("time (ms)")
ylabel("P")

yline(.05, '--','p=.05','Color','black','HandleVisibility','off')
legend
set(gcf,'color','white')
fontsize(30,'points')
hold off

%%
figure;
hold on
plot(time_lags, p_Interaction_envelope, 'DisplayName','Interaction Pval','LineWidth',2.5,'Color',[.5 .8 .5])
plot(time_lags, p_Performance_envelope, 'DisplayName','Performance Condition Pval','LineWidth',2.5,'Color',[.8 0 .5])
plot(time_lags, p_Tempo_envelope, 'DisplayName','Tempo Condition Pval','LineWidth',2.5,'Color',[1 .6 0])

title("Envelope Pvals across Time Lags")
xlabel("time (ms)")
ylabel("P")

yline(.05, '--','p=.05','Color','black','HandleVisibility','off')
set(gcf,'color','white')
fontsize(30,'points')
legend
hold off


%% Scatter Plot with questionairre
questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction (Reverse-Scored)','Spontaneity'};
times = 14:26;
electrodes = 1:29;
figure;
t = tiledlayout('flow');
t.Title.String = "Questionaire vs TRF Corr (avg across 100-200 ms and Electrodes)";
t.Title.FontWeight = 'Bold';
for Q_num = 1:8
    nexttile;
    q_TRF_flux = squeeze(mean(all_flux_TRFs(:,:,electrodes,times),[3,4]));
    Q_data = squeeze(aggregate(:,:,Q_num));
    
    scatter(Q_data(1:end),q_TRF_flux(1:end));
    title(questions(Q_num))
    xlabel("Rating 1-5")
    ylabel("TRF Corr")
    xlim([1, 5])
    ylim([-.03, .05])
    lsline;
    r = corrcoef(Q_data(1:end)',q_TRF_flux(1:end)');
    r_squared = r(2,1)^2;
    text(4, .04,"r^2: "+num2str(round(r_squared,3)))
    if r_squared > .1
        error("success?")
    end
    set(gcf,'color','white')
    fontsize(20,'points')
end

%% Plot flux TRFs by Q-score
questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction (Reverse-Scored)','Spontaneity'};
conditions = 1:4;
electrodes = [4, 8, 26];
figure;
t = tiledlayout('flow');
t.Title.String = "TRF Corr by Questionnaire High/Low Rating";
% t.Title.FontWeight = 'Bold';
for Q_num = 1:8
    if Q_num == 3
        nexttile;
        Q_data = squeeze(aggregate(:,conditions,Q_num));
        q_low = find(Q_data<mean(Q_data,'all'));
        q_high = find(Q_data>mean(Q_data,'all')); % mean split
        n_low = length(q_low);
        n_high = length(q_high);
        [rlow, clow] = ind2sub([21,4],q_low);
        [rhigh, chigh] = ind2sub([21,4],q_high);
        low_TRF_flux = zeros(length(q_low),size(all_flux_TRFs,4));
        high_TRF_flux = zeros(length(q_high),size(all_flux_TRFs,4));
        numLow = 0;
        numHigh = 0;
        for i = 1:length(q_low)
            % if (clow(i) == 3 || clow(i) == 4) %% change the conditions (1,2) for live or (3,4) for rec
                numLow = numLow + 1;
                low_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rlow(i),clow(i),electrodes,:),3));
            % end
        end
        for i = 1:length(q_high)
            % if (chigh(i) == 3 || chigh(i) == 4)
                numHigh = numHigh + 1;
                high_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rhigh(i),chigh(i),electrodes,:),3));
            % end
        end
        low_TRF_flux = low_TRF_flux(all(low_TRF_flux,2),:);
        high_TRF_flux = high_TRF_flux(all(high_TRF_flux,2),:);
        errlow = std(low_TRF_flux,0,1)/sqrt(length(low_TRF_flux));
        errhigh = std(high_TRF_flux,0,1)/sqrt(length(high_TRF_flux));
        
        hold on
        shadedErrorBar(time_lags,mean(high_TRF_flux,1),errhigh,'lineprops',{'Color',[.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(numHigh)})
        shadedErrorBar(time_lags,mean(low_TRF_flux,1),errlow,'lineprops',{'Color',[0, .2, .8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(numLow)})
        yline(0,'--','LineWidth',3,'Color','black','LineWidth', 1,'HandleVisibility','off')
        xline(164.0625,'label','164.06 ms','LineWidth', 1,'color',[0.8, 0.1, 0.2],'LabelVerticalAlignment', 'bottom','HandleVisibility','off') 
        hold off;
        title(questions(Q_num))
        xlabel("Time Lag (ms)")
        ylabel("TRF Corr")
        ylim([-.03, .05])
        legend;
    end
end

nexttile;
Q_data = squeeze(mean(aggregate(:,conditions,:),3));
q_low = find(Q_data<mean(Q_data,'all'));
q_high = find(Q_data>mean(Q_data,'all'));
n_low = length(q_low);
n_high = length(q_high);
[rlow, clow] = ind2sub([21,4],q_low);
[rhigh, chigh] = ind2sub([21,4],q_high);
low_TRF_flux = zeros(length(q_low),size(all_flux_TRFs,4));
high_TRF_flux = zeros(length(q_high),size(all_flux_TRFs,4));
numLow = 0;
numHigh = 0;
for i = 1:length(q_low)
    % if (clow(i) == 3 || clow(i) == 4)
        numLow = numLow + 1;
        low_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rlow(i),clow(i),electrodes,:),3));
    % end
end
for i = 1:length(q_high)
    % if (chigh(i) == 3 || chigh(i) == 4)
        numHigh = numHigh + 1;
        high_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rhigh(i),chigh(i),electrodes,:),3));
    % end
end
low_TRF_flux = low_TRF_flux(all(low_TRF_flux,2),:);
high_TRF_flux = high_TRF_flux(all(high_TRF_flux,2),:);
% 
errlow = std(low_TRF_flux,0,1)/sqrt(length(low_TRF_flux));
errhigh = std(high_TRF_flux,0,1)/sqrt(length(high_TRF_flux));

hold on
shadedErrorBar(time_lags,mean(high_TRF_flux),errhigh,'lineprops',{'Color',[0.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(numHigh)})
shadedErrorBar(time_lags,mean(low_TRF_flux),errlow,'lineprops',{'Color',[0, .2, 0.8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(numLow)})
yline(0,'--','LineWidth',3,'Color','black','LineWidth', 1,'HandleVisibility','off')
xline(164.0625,'label','164.06 ms','LineWidth', 1,'color',[0.8, 0.1, 0.2],'LabelVerticalAlignment', 'bottom','HandleVisibility','off') 
hold off;
title("Aggregate")
xlabel("Time Lag (ms)")
ylabel("TRF Corr")
ylim([-.03, .05])
legend;
fontsize(20,'points')
set(gcf,'color','white')
%% Unequal sample t test for above plot
time_lag_index = 22;
[h, p] = ttest2(mean(low_TRF_flux(:, time_lag_index), 2), mean(high_TRF_flux(:, time_lag_index), 2), 'Vartype', 'unequal');



%% Plot envelope TRFs by Q-score
questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction (Reverse-Scored)','Spontaneity'};
electrodes = 1:29;
figure;
t = tiledlayout('flow');
t.Title.String = "TRF Corr by Questionnaire High/Low Rating";
t.Title.FontWeight = 'Bold';
for Q_num = 1:8
    nexttile;
    Q_data = squeeze(aggregate(:,:,Q_num));
    q_low = find(Q_data<=3);
    q_high = find(Q_data>=3);
    n_low = length(q_low);
    n_high = length(q_high);
    [rlow, clow] = ind2sub([21,4],q_low);
    [rhigh, chigh] = ind2sub([21,4],q_high);
    low_TRF_envelope = zeros(length(q_low),size(all_envelope_TRFs,4));
    high_TRF_envelope = zeros(length(q_high),size(all_envelope_TRFs,4));
    for i = 1:length(q_low)
        low_TRF_envelope(i,:) = squeeze(mean(all_envelope_TRFs(rlow(i),clow(i),electrodes,:),3));
    end
    for i = 1:length(q_high)
        high_TRF_envelope(i,:) = squeeze(mean(all_envelope_TRFs(rhigh(i),chigh(i),electrodes,:),3));
    end

    errlow = std(low_TRF_envelope,0,1)/sqrt(length(low_TRF_envelope));
    errhigh = std(high_TRF_envelope,0,1)/sqrt(length(high_TRF_envelope));
    
    hold on
    shadedErrorBar(time_lags,mean(high_TRF_envelope,1),errhigh,'lineprops',{'Color',[.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(n_high)})
    shadedErrorBar(time_lags,mean(low_TRF_envelope,1),errlow,'lineprops',{'Color',[0, .2, .8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(n_low)})

    hold off;
    title(questions(Q_num))
    xlabel("time(ms)")
    ylabel("TRF Corr")
    ylim([-.03, .05])
    legend;
    fontsize(15,'points')
end

nexttile;
Q_data = squeeze(mean(aggregate,3));
q_low = find(Q_data<3);
q_high = find(Q_data>3);
n_low = length(q_low);
n_high = length(q_high);
[rlow, clow] = ind2sub([21,4],q_low);
[rhigh, chigh] = ind2sub([21,4],q_high);
low_TRF_envelope = zeros(length(q_low),size(all_envelope_TRFs,4));
high_TRF_envelope = zeros(length(q_high),size(all_envelope_TRFs,4));
for i = 1:length(q_low)
    low_TRF_envelope(i,:) = squeeze(mean(all_envelope_TRFs(rlow(i),clow(i),electrodes,:),3));
end
for i = 1:length(q_high)
    high_TRF_envelope(i,:) = squeeze(mean(all_envelope_TRFs(rhigh(i),chigh(i),electrodes,:),3));
end

% 
errlow = std(low_TRF_envelope,0,1)/sqrt(length(low_TRF_envelope));
errhigh = std(high_TRF_envelope,0,1)/sqrt(length(high_TRF_envelope));

hold on
shadedErrorBar(time_lags,mean(high_TRF_envelope),errhigh,'lineprops',{'Color',[0.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(n_high)})
shadedErrorBar(time_lags,mean(low_TRF_envelope),errlow,'lineprops',{'Color',[0, .2, 0.8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(n_low)})
hold off;
title("Aggregate")
xlabel("time(ms)")
ylabel("TRF Corr")
ylim([-.03, .05])
legend;
fontsize(15,'points')
set(gcf,'color','white')
%% PLot by condition and Q score
conditions = 1:4;
electrodes = [4, 8, 26];
figure;
% t = tiledlayout('flow');
% t.Title.String = "PLV by Questionnaire High/Low Rating";
% t.Title.FontWeight = 'Bold';


Q_data = squeeze(mean(aggregate(:,conditions,:),3));
q_low = find(Q_data<median(Q_data(:, [1 3]),'all'));
q_high = find(Q_data>median(Q_data(:, [1 3]),'all'));
n_low = length(q_low);
n_high = length(q_high);
[rlow, clow] = ind2sub([21,4],q_low);
[rhigh, chigh] = ind2sub([21,4],q_high);
low_flux_cond = zeros(length(q_low),size(all_flux_TRFs,4));
high_flux_cond = zeros(length(q_high),size(all_flux_TRFs,4));
numLow = 0;
numHigh = 0;
for i = 1:length(q_low)
    if ((clow(i) == 1) || (clow(i)==3)) % check condition is live or rec
        numLow = numLow + 1;
        low_flux_cond(i,:) = squeeze(mean(all_flux_TRFs(rlow(i),clow(i),electrodes,:),3));
    end
end
for i = 1:length(q_high)
    if ((chigh(i) == 1) || (chigh(i)==3))
        numHigh = numHigh + 1;
        high_flux_cond(i,:) = squeeze(mean(all_flux_TRFs(rhigh(i),chigh(i),electrodes,:),3));
    end
end
% test_PLV = low_PLV(all(low_PLV,2),:);
% high_TRF_flux = high_TRF_flux(all(high_TRF_flux,2),:);
% 
errlow = std(low_flux_cond,0,1)/sqrt(length(low_flux_cond));
errhigh = std(high_flux_cond,0,1)/sqrt(length(high_flux_cond));

hold on
shadedErrorBar(time_lags,mean(high_flux_cond),errhigh,'lineprops',{'Color',[0.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(numHigh)})
shadedErrorBar(time_lags,mean(low_flux_cond),errlow,'lineprops',{'Color',[0, .2, 0.8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(numLow)})
hold off;
title("Aggregate")
xlabel("Frequency [Hz]")
ylabel("TRF Corr")
ylim([-.03, .05])
legend;
fontsize(15,'points')
set(gcf,'color','white')
%% plot flux by diff Q score Live/Rec see PLV_mirCochlear --> inconclusive results
figure;
t = tiledlayout('flow');
t.Title.String = 'PLVs by Condition and Q-diff';
% hold on

for condition = [1 3 2 4]
    switch condition
        case 1
            % condition = [1 2]
            nexttile; hold on;
            colorHigh = [0.8, .2, 0];
            colorLow = [0.6, .2, 0];
            label = "Live Fast";
            title("Fast")
            % xlim([0.2,20.2])
            % ylim([0, .14])
            % yticks(0:.02:.14)
            % xticks(sort([xticks 8.4])); xline([8.4],"HandleVisibility","off","LineWidth",1);
            % xlabel("Frequency [Hz]"); ylabel("PLV");
        case 3
            colorHigh = [0, .2, 0.8];
            colorLow = [0, .2, 0.4];
            label = "Recorded Fast";
        case 2
            % condition = [3 4]
            hold off; nexttile; hold on;
            colorHigh = [0.8, .2, 0];
            colorLow = [0.6, .2, 0];
            label = "Live Slow";
            title("Slow")
            % xlim([0.2,6.2])
            % ylim([0, .14])
            % yticks(0:.02:.14)
            % xticks(sort([xticks .8])); xline([.8],"HandleVisibility","off","LineWidth",1);
            % xlabel("Frequency [Hz]"); ylabel("PLV");
    
        case 4
            colorHigh = [0, .2, 0.8];
            colorLow = [0, .2, 0.4];
            label = "Recorded Slow";
    end
    disp(condition)
    lowflux = squeeze(mean(all_flux_TRFs(participant_ranking(1:10,1),condition,electrodes,:),[2 3]));
    highflux = squeeze(mean(all_flux_TRFs(participant_ranking(11:21,1),condition,electrodes,:),[2 3]));
    
    errlow = 2*std(lowflux,0,1)/sqrt(length(lowflux));
    errhigh = 2*std(highflux,0,1)/sqrt(length(highflux));
    
    shadedErrorBar(time_lags,mean(highflux),errhigh,'lineprops',{'Color',colorHigh,'LineWidth',2,'DisplayName',label + " High Rating Diff"})
    shadedErrorBar(time_lags,mean(lowflux),errlow,'lineprops',{'Color',colorLow,'LineWidth',2,'DisplayName',label + " Low Rating Diff"})
    legend;
end


hold off;
xlabel("Frequency [Hz]")
ylabel("PLV")

% legend;
fontsize(20,'points')
set(gcf,'color','white')


%% Plot TERTILE flux TRFs by Aggregate Q-score 
questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction (Reverse-Scored)','Spontaneity'};
groupings = {'LiveFast', 'LiveSlow','RecFast','RecSlow'};
all_conditions = {1,2,3,4};
electrodes = [4, 8, 26];
figure;
t = tiledlayout('flow');
t.Title.String = "TRF Corr by Aggregate Questionnaire Score (Tertile Split)";
t.Title.FontWeight = 'Bold';
for group = 1:4
    conditions = all_conditions(group);
    nexttile;
    Q_data = squeeze(mean(aggregate(:,conditions{1},:),3));
    thresholds = quantile(reshape(Q_data,numel(Q_data),1),[0.333, 0.666]);
    q_low = find(Q_data<thresholds(1));
    q_med = find(Q_data>thresholds(1) & Q_data<thresholds(2));
    q_high = find(Q_data>thresholds(2)); % tertile split
    n_low = length(q_low);
    n_med = length(q_med);
    n_high = length(q_high);
    [rlow, clow] = ind2sub([21,4],q_low);
    [rmed, cmed] = ind2sub([21,4],q_med);
    [rhigh, chigh] = ind2sub([21,4],q_high);
    low_TRF_flux = zeros(length(q_low),size(all_flux_TRFs,4));
    med_TRF_flux = zeros(length(q_med),size(all_flux_TRFs,4));
    high_TRF_flux = zeros(length(q_high),size(all_flux_TRFs,4));
    for i = 1:length(q_low)
        low_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rlow(i),clow(i),electrodes,:),3));
    end
    for i = 1:length(q_med)
        med_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rmed(i),cmed(i),electrodes,:),3));
    end
    for i = 1:length(q_high)
        high_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rhigh(i),chigh(i),electrodes,:),3));
    end
    
    % 
    errlow = std(low_TRF_flux,0,1)/sqrt(length(low_TRF_flux));
    errmed = std(med_TRF_flux,0,1)/sqrt(length(med_TRF_flux));
    errhigh = std(high_TRF_flux,0,1)/sqrt(length(high_TRF_flux));
    hold on
    shadedErrorBar(time_lags,mean(high_TRF_flux,1),errhigh,'lineprops',{'Color',[.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(n_high)})
    shadedErrorBar(time_lags,mean(med_TRF_flux,1),errmed,'lineprops',{'Color',[.1, .8, .1],'LineWidth',2,'DisplayName',"Medium Rating n=" + num2str(n_med)})
    shadedErrorBar(time_lags,mean(low_TRF_flux,1),errlow,'lineprops',{'Color',[0, .2, .8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(n_low)})
    hold off;
    title(groupings(group))
    xlabel("time lag (ms)")
    ylabel("TRF Corr")
    ylim([-.03, .05])
    legend;
    fontsize(15,'points')
    set(gcf,'color','white')
end
%% Plot Median flux TRFs by Aggregate Q-score 
questions = {'Engagement','Enjoyment','Familiarity','Pleasure','Focus','Investment','Distraction (Reverse-Scored)','Spontaneity'};
groupings = {'LiveFast', 'LiveSlow','RecFast','RecSlow'};
all_conditions = {1,2,3,4};
electrodes = [4, 8, 26];
figure;
t = tiledlayout('flow');
t.Title.String = "TRF Corr by Aggregate Questionnaire Score (Median Split)";
t.Title.FontWeight = 'Bold';
for group = 1:4
    conditions = all_conditions(group);
    nexttile;
    Q_data = squeeze(mean(aggregate(:,conditions{1},:),3));
    threshold = median(Q_data,"all");
    q_low = find(Q_data<thresholds);
    q_high = find(Q_data>threshold); % median split
    n_low = length(q_low);
    n_high = length(q_high);
    [rlow, clow] = ind2sub([21,4],q_low);
    [rhigh, chigh] = ind2sub([21,4],q_high);
    low_TRF_flux = zeros(length(q_low),size(all_flux_TRFs,4));
    high_TRF_flux = zeros(length(q_high),size(all_flux_TRFs,4));
    for i = 1:length(q_low)
        low_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rlow(i),clow(i),electrodes,:),3));
    end
    for i = 1:length(q_high)
        high_TRF_flux(i,:) = squeeze(mean(all_flux_TRFs(rhigh(i),chigh(i),electrodes,:),3));
    end
    
    % 
    errlow = std(low_TRF_flux,0,1)/sqrt(length(low_TRF_flux));
    errhigh = std(high_TRF_flux,0,1)/sqrt(length(high_TRF_flux));
    hold on
    shadedErrorBar(time_lags,mean(high_TRF_flux,1),errhigh,'lineprops',{'Color',[.8, .2, 0],'LineWidth',2,'DisplayName',"High Rating n=" + num2str(n_high)})
    shadedErrorBar(time_lags,mean(low_TRF_flux,1),errlow,'lineprops',{'Color',[0, .2, .8],'LineWidth',2,'DisplayName',"Low Rating n=" + num2str(n_low)})
    hold off;
    title(groupings(group))
    xlabel("time lag (ms)")
    ylabel("TRF Corr")
    ylim([-.03, .05])
    legend;
    fontsize(15,'points')
    set(gcf,'color','white')
end





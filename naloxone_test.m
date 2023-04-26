%% Specify signature to test and load it
test_reward = false; % default: false
load([resultsdir filesep 'pos_affect_pattern.mat']); %initialize with proper voxels included

if test_reward
    pos_affect_pattern = brs;
end

%% Load brain responses (BOLD) to erotic images

files = dir([datadir filesep 'MRI' filesep 'fMRI' filesep 'Saline_greater_than_Naloxone' filesep 'erotic' filesep 'source' filesep '*.nii']);
P={};
for f=1:length(files)
    P{f}=[files(f).folder filesep files(f).name]; 

end

saline_vs_naloxone_erotic = fmri_data(P); %load data
outliers_erotic = plot(saline_vs_naloxone_erotic); %look for potential outliers

%% Compute pattern expression and use the bootstrap for stats on responses to erotic images

pexp_erotic = apply_mask(saline_vs_naloxone_erotic,pos_affect_pattern,'pattern_expression','cosine_similarity');
bs_erotic = bootstrp(10000,@mean,pexp_erotic(~outliers_erotic));
bs_Z_erotic = mean(bs_erotic)./std(bs_erotic);
b_P_erotic = 2*normcdf(-1*abs(bs_Z_erotic),0,1);
bs_erotic_es = bootstrp(10000,@boot_eff,pexp_erotic);
%% Load brain responses (BOLD) to monetary rewards

files = dir([datadir filesep 'MRI' filesep 'fMRI' filesep 'Saline_greater_than_Naloxone' filesep 'money' filesep 'source' filesep '*.nii']);
P={};
for f=1:length(files)
    P{f}=[files(f).folder filesep files(f).name];

end

saline_vs_naloxone_monetary = fmri_data(P); %load data
outliers_money = plot(saline_vs_naloxone_monetary); %look for potential outliers

%% Compute pattern expression and use the bootstrap for stats on responses to monetary rewards

pexp_money = apply_mask(saline_vs_naloxone_monetary,pos_affect_pattern,'pattern_expression','cosine_similarity');
bs_monetary = bootstrp(10000,@mean,pexp_money(~outliers_money));
bs_monetary_es = bootstrp(10000,@boot_eff,pexp_money);
b_P_monetary = 2*normcdf(-1*abs(mean(bs_monetary)./std(bs_monetary)),0,1);

%% 

bs_means = bootstrp(10000,@mean,[pexp_erotic(~outliers_money & ~outliers_erotic) pexp_money(~outliers_money & ~outliers_erotic)]);
figure;
distributionPlot([bs_means(:,1),bs_means(:,2)],'showMM',5,'color',[.6 .6 .6]);
hold on;axis square;
hl = findobj(gca,'type','line');
set(hl,'Linewidth',2,'color','k')
ylabel '\Delta Pattern Expression (Placebo vs. Naloxone)'
set(gca,'XTickLabel',{'Erotic Images','Money'})
drawnow;snapnow;


%% Load behavioral data and conduct the same statistical tests
load([datadir filesep 'Ratings' filesep 'ratings.mat')
behav_pleasantness_erotic = [(behav([5]).vals(1:19)+behav([7]).vals(1:19))/2 (behav([5]).vals(20:end)+behav([7]).vals(20:end))/2];
behav_pleasantness_money = [(behav([1]).vals(1:19)+behav([3]).vals(1:19))/2 (behav([1]).vals(20:end)+behav([3]).vals(20:end))/2];

bs_erotic_behav = bootstrp(10000,@mean,behav_pleasantness_erotic(:,2)-behav_pleasantness_erotic(:,1));
b_P_erotic_behav = 2*normcdf(-1*abs(mean(bs_erotic_behav)./std(bs_erotic_behav)),0,1);

bs_money_behav = bootstrp(10000,@mean,behav_pleasantness_money(:,2)-behav_pleasantness_money(:,1));
b_P_money_behav = 2*normcdf(-1*abs(mean(bs_money_behav)./std(bs_money_behav)),0,1);

figure;
bs_means_behav = bootstrp(10000,@mean,[behav_pleasantness_erotic(:,2)-behav_pleasantness_erotic(:,1) (behav_pleasantness_money(:,2)-behav_pleasantness_money(:,1))]);
distributionPlot([bs_means_behav(:,1),bs_means_behav(:,2)],'showMM',5,'color',[.6 .6 .6]);
hl = findobj(gca,'type','line');
set(hl,'Linewidth',2,'color','k')

hold on;axis square;

ylabel '\Delta Pleasure Rating (Placebo vs. Naloxone)'
set(gca,'XTickLabel',{'Erotic Images','Money'})
drawnow;snapnow;

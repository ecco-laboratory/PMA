%%
load('pos_affect_pattern.mat')

files = dir('C:\Users\phili\OneDrive - Emory University\Projects\PositiveValenceMegaAnalysis\Data\ReceptorABA\*.gz');
clear P
for f=1:length(files); 
    P{f}=[files(f).folder filesep files(f).name];
end

nt_maps = fmri_data(P);
nt_maps = preprocess(nt_maps,'smooth',6);


leg = {'DRD1','DRD2','DRD3','OPRD1','OPRK1','OPRM1'};
for ii=1:6; 
    tv = nt_maps; tv.dat=tv.dat(:,ii); tv.fullpath = [datadir filesep 'neurotransmitter_maps_' leg{ii} '.nii'];
% write(tv,'overwrite');
end

%%
load('pos_affect_pattern.mat');
load bs_b
tv = nt_maps;
tv = apply_mask(tv,pos_affect_pattern);
tv = resample_space(tv,pos_affect_pattern);
tv = replace_empty(tv);
%%
for rgr=2:19

for k=1:size(bs_b,1)


load('pos_affect_pattern.mat');


pos_affect_pattern.dat = squeeze(bs_b(k,2:end,rgr))';
pos_affect_pattern=replace_empty(pos_affect_pattern);


%spatial regression
[b(k,:), dev, st]=glmfit(tv.dat(tv.dat(:,1)~=0 & pos_affect_pattern.dat~=0,:),pos_affect_pattern.dat(tv.dat(:,1)~=0 & pos_affect_pattern.dat~=0));

R = corr((b(k,2:end)*(tv.dat(tv.dat(:,1)~=0 & pos_affect_pattern.dat~=0,:))')'+b(k,1),pos_affect_pattern.dat(tv.dat(:,1)~=0 & pos_affect_pattern.dat~=0));
for i=1:6
    rs(k,rgr,i)= corr((b(k,2:end)*(tv.dat(tv.dat(:,1)~=0 & pos_affect_pattern.dat~=0,:))')'+b(k,1),(tv.dat(tv.dat(:,i)~=0 & pos_affect_pattern.dat~=0,i)));
end

stats = image_similarity_plot(pos_affect_pattern, 'correlation', 'mapset',tv,'notable','nofigure','noplot');
r_nt_bootstrap(k,rgr,:)=stats.r;
k/3000
end
end


for i=1:6
    for j=1:6
            mean_bs = mean(b(:,1+i)-b(:,1+j));
            std_bs = std(b(:,1+i)-b(:,1+j));
            bs_Z_paired(i,j)=mean_bs./std_bs;
            bs_P_paired(i,j) = 2*normcdf(-1*abs(bs_Z_paired(i,j)),0,1);
    end
end


mean_bs = mean(b(:,2:end));
std_bs = std(b(:,2:end));
bs_Z=mean_bs./std_bs;
bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    


mean_bs_diff = mean(b(:,end)-mean(b(:,2:end-1),2));
std_bs_diff = std(b(:,end)-mean(b(:,2:end-1),2));
bs_Z_diff=mean_bs_diff./std_bs_diff;
bs_P_diff = 2*normcdf(-1*abs(bs_Z_diff),0,1);

%% plot regression
C = seaborn_colors(10);
figure;

distributionPlot(b(:,2:end),'showMM',0,'color',C([4 3 6 1 5 2],:)); legend(leg);
set(gca,'XTick','')
ylabel 'Beta Coefficient'

%% correlations
figure;

distributionPlot(squeeze(r_nt_bootstrap(:,2,:)),'showMM',0,'color',C([4 3 6 1 5 2],:)); legend(leg);
set(gca,'XTick','')
ylabel 'Correlation Coefficient'

mean_bs = mean(r_nt_bootstrap);
std_bs = std(r_nt_bootstrap);
bs_Z=mean_bs./std_bs;
bs_P = 2*normcdf(-1*abs(bs_Z),0,1)

%%
figure;imagesc(squeeze(mean(r_nt_bootstrap(:,2:end,:))))
colorbar
set(gca,'XTickLabel',leg)
set(gca,'YTick',1:18);
set(gca,'YTickLabel',{'Positive Affect','Pain','Cognitive Control','Negative Affect','Music','Food','Erotic','Monetary','Social','Thermal','Visceral','Mechanical','Working Memory','Response Selection','Response Conflict','Visual','Social','Auditory'})
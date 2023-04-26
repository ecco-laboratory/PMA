% Plot results of searchlight analyis and their correspondence with Buckner lab parcellation


%% Make riverplot of searchlight maps and Buckner lab resting-state parcellation
load([resultsdir filesep 'pos_affect_pattern.mat'])
load([tempdir filesep 'all_regs_by_roi.mat'])

names = {'vmPFC'  'dmPFC' 'sgACC' 'pMCC' 'pACC' 'aMCC' subcort.labels{:} insula.labels{:}};

tv = condf2indic(all_regs_by_roi.dat);

all_regs_by_roi.dat = tv(:,2:7);
image_obj3= resample_space(all_regs_by_roi,pos_affect_pattern,'nearest');
layer1colors = {.8*[1 .5 0] [0 .2 1]};
% C={[120 18 134]/255 [70 130 180]/255  [0 118 14]/255  [196 58 250]/255  [220 248 164]/255  [230 148 34]/255  [205 62 78]/255 };

newdat = zeros(size(pos_affect_pattern.dat,1),2);
newdat(:,1) = pos_affect_pattern.dat.*double(pos_affect_pattern.dat>0);
newdat(:,2) = -1*pos_affect_pattern.dat.*double(pos_affect_pattern.dat<0);
pos_affect_pattern.dat=newdat;

figure;
image_obj3.image_names=names(1:6);
layer2colors = seaborn_colors(length(image_obj3.image_names));

pos_affect_pattern.image_names = {'Positive Coefficients' 'Negative Coefficients'};
riverplot(pos_affect_pattern,'layer2',image_obj3,'pos','thin','colors1',layer1colors,'colors2',layer2colors,'reorder')
drawnow;snapnow;

%%
figure;
all_regs_by_roi.dat = tv(:,[8:10 12:23]);
image_obj3= resample_space(all_regs_by_roi,pos_affect_pattern,'nearest');
image_obj3.image_names=names([8:10 12:23]-1);
layer2colors = seaborn_colors(length(image_obj3.image_names));

pos_affect_pattern.image_names = {'Positive Coefficients' 'Negative Coefficients'};
riverplot(pos_affect_pattern,'layer2',image_obj3,'pos','thin','colors1',layer1colors,'colors2',layer2colors,'reorder')
drawnow;snapnow;

%%
all_regs_by_roi.dat = tv(:,[11 24:27]);
image_obj3= resample_space(all_regs_by_roi,pos_affect_pattern,'nearest');
image_obj3.image_names=names([11 24:27]-1);
layer2colors = seaborn_colors(length(image_obj3.image_names));

pos_affect_pattern.image_names = {'Positive Coefficients' 'Negative Coefficients'};
riverplot(pos_affect_pattern,'layer2',image_obj3,'pos','thin','colors1',layer1colors,'colors2',layer2colors,'reorder')
drawnow;snapnow;


%%
all_regs_by_roi.dat = tv(:,28:end);
all_regs_by_roi.dat = all_regs_by_roi.dat(:,1:2:end)+all_regs_by_roi.dat(:,2:2:end);

image_obj3= resample_space(all_regs_by_roi,pos_affect_pattern,'nearest');

image_obj3.image_names=names(27:2:end);
image_obj3.image_names = strrep(image_obj3.image_names,'L','');
image_obj3.image_names = strrep(image_obj3.image_names,'Ctx','');


layer2colors = seaborn_colors(length(image_obj3.image_names));

pos_affect_pattern.image_names = {'Positive Coefficients' 'Negative Coefficients'};
riverplot(pos_affect_pattern,'layer2',image_obj3,'pos','thin','colors1',layer1colors,'colors2',layer2colors,'reorder')
drawnow;snapnow;

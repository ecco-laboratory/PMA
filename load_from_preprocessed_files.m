%% MUSIC 

%kragel & labar 2015 amusing music
load([datadir filesep 'Kragel2015EmotionClassification.mat'])
load([datadir filesep 'nws.mat'])
dat=remove_empty(dat);
dat.Y(dat.removed_images)=[];

tmusic_dat=dat;
iM=iM(7:end);
tmusic_dat.Y=tmusic_dat.Y(iM==1);
tmusic_dat.dat=tmusic_dat.dat(:,iM==1);
subind = S(IB);
subind=subind(7:end);
subind=subind(iM==1);
music_dat=tmusic_dat;


for s=1:32
music_dat.dat(:,s)=nanmean(tmusic_dat.dat(:,tmusic_dat.Y==2 & subind'==s),2);
end

music_dat.dat=music_dat.dat(:,1:32);

for i=1:size(music_dat.dat,2) 
    music_dat.dat(isnan(music_dat.dat(:,i)),i)=nanmean(music_dat.dat(:,i));
end
music_dat.removed_images=0;


% pleasant music: 171 con - 2; kragel & labar 2015;
%files = dir('C:\Users\phili\ds000171-download\*control*\**\con_0001.nii');
%for f=1:length(files)
%    P{f} = [files(f).folder filesep files(f).name];
%end
pleasant_music = fmri_data([datadir filesep 'Study2.nii']);


%% FOOD 
% pleasant food images: 2270 - con 1; 157 - con 1; 
files = dir('C:\Users\phili\ds002270-download\*sub*\**\con_0001.nii');
clear P;
for f=1:length(files)
    P{f} = [files(f).folder filesep files(f).name];
end
pleasant_food = fmri_data(P);

% pleasant food images: 2270 - con 1; 157 - con 1; 
files = dir('C:\Users\phili\ds000157-download\*sub*\**\con_0001.nii');
clear P;
for f=1:length(files)
    P{f} = [files(f).folder filesep files(f).name];
end
pleasant_food_two = fmri_data(P);

%% SEX 
% positive iaps images/sexy images:  1491 - con 1; kragel, reddan et al. 2019
sexual_iaps = fmri_data('C:\Users\phili\Dropbox (Cognitive and Affective Neuroscience Laboratory)\Cooperation Valencespace\data\neuroimaging data\IAPS_firstlevel_contrasts_sexual_positive_negative\IAPS_n18_sexual.nii');

files = dir('C:\Users\phili\ds001491-download\*sub*\**\con_0001.nii');
clear P;
for f=1:length(files)
    P{f} = [files(f).folder filesep files(f).name];
end
sexual_images = fmri_data(P);

%% MONEY

files = dir('C:\Users\phili\ds002041-download\**\con_0001.nii');
clear P;
for f=1:length(files)
    P{f} = [files(f).folder filesep files(f).name];
end

td_reward_cue =  fmri_data(P);


files = dir('C:\Users\phili\ds000005-download\**\con_0001.nii');
clear P;
for f=1:length(files)
    P{f} = [files(f).folder filesep files(f).name];
end

gain_cue = fmri_data(P);

%% SOCIAL

% neonate videos
own_baby = fmri_data('C:\Users\phili\ds003136-download\own_positive.nii');

% social cues (Tye/Saxe)
files = dir('C:\Users\phili\ds003242-download\**\con_0001.nii');
clear P;
for f=1:length(files)
    P{f} = [files(f).folder filesep files(f).name];
end
social = fmri_data(P);

%% concatenate

music_dat = resample_space(music_dat,pleasant_food);
pleasant_music = resample_space(pleasant_music,pleasant_food);

pleasant_food_two = resample_space(pleasant_food_two,pleasant_music);

sexual_images = resample_space(sexual_images,pleasant_music);
sexual_iaps = resample_space(sexual_iaps,pleasant_music);

td_reward_cue = resample_space(td_reward_cue,pleasant_music);
gain_cue = resample_space(gain_cue,pleasant_music);

own_baby = resample_space(own_baby,pleasant_music);
social  = resample_space(social,pleasant_music);

alldat = music_dat;
alldat.dat = [alldat.dat pleasant_music.dat pleasant_food.dat pleasant_food_two.dat sexual_images.dat sexual_iaps.dat td_reward_cue.dat gain_cue.dat own_baby.dat social.dat];
study = [ones(size(music_dat.dat,2),1); 2*ones(size(pleasant_music.dat,2),1); 3*ones(size(pleasant_food.dat,2),1); 4*ones(size(pleasant_food_two.dat,2),1); 5*ones(size(sexual_images.dat,2),1); 6*ones(size(sexual_iaps.dat,2),1); 7*ones(size(td_reward_cue.dat,2),1); 8*ones(size(gain_cue.dat,2),1) ; 9*ones(size(own_baby.dat,2),1) ; 10*ones(size(social.dat,2),1)];
    

%% add pain, cognitive control, and negative emotion

data_obj = load_image_set('kragel18_alldata');
nstudy = max(study);
for i=1:18
    study = [study; nstudy+i*ones(15,1)]; 
end
data_obj=resample_space(data_obj,alldat);
alldat.dat = [alldat.dat data_obj.dat];

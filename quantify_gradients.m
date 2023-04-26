nac = select_atlas_subset(load_atlas('CIT168'),{'NAC'});
VeP = select_atlas_subset(load_atlas('CIT168'),{'VeP'});
insula = select_atlas_subset(load_atlas('glasser'),{'Ctx_PoI1_L','Ctx_Ig_L','Ctx_PI_L','Ctx_POI2_L','Ctx_MI_L','Ctx_AI_L','Ctx_AVI_L','Ctx_AAIC_L','Ctx_PoI1_R','Ctx_Ig_R','Ctx_PI_R','Ctx_POI2_R','Ctx_MI_R','Ctx_AI_R','Ctx_AVI_R','Ctx_AAIC_R'});
load([tempdir filesep 'bs_b.mat'])
load([resultsdir filesep 'pos_affect_pattern.mat'])
tv = pos_affect_pattern;
%%


for i=1:3000
    pos_affect_pattern = tv;
    pos_affect_pattern.dat = bs_b(i,2:end,2)';
    pos_affect_pattern = apply_mask(pos_affect_pattern,nac);


    b_nac(i,:) = glmfit(abs(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,1)-41),pos_affect_pattern.dat);


    pos_affect_pattern = tv;
    pos_affect_pattern.dat = bs_b(i,2:end,2)';
    pos_affect_pattern = apply_mask(pos_affect_pattern,VeP);
    b_VeP(i,:) = glmfit(abs(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,2)),pos_affect_pattern.dat);

    pos_affect_pattern = tv;
    pos_affect_pattern.dat = bs_b(i,2:end,2)';
    pos_affect_pattern = apply_mask(pos_affect_pattern,insula);
    b_insula(i,:) = glmfit(abs(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,2)),pos_affect_pattern.dat);

    st=tpaps([abs(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,1)-41),(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,2))]',pos_affect_pattern.dat');
    points = fnplt(st);
    [u_insula(i,:,:),v_insula(i,:,:)] = gradient(-1*downsample((downsample(points{3},5))',5));

end



mean_b_nac = mean(1000*b_nac); %rescale arbitrary units
std_b_nac = std(1000*b_nac);
b_Z_nac=mean_b_nac./std_b_nac;
b_P_nac = 2*normcdf(-1*abs(b_Z_nac),0,1);



mean_b_VeP = mean(1000*b_VeP);
std_b_VeP = std(1000*b_VeP);
b_Z_VeP=mean_b_VeP./std_b_VeP;
b_P_VeP = 2*normcdf(-1*abs(b_Z_VeP),0,1);



mean_b_insula = mean(1000*b_insula);
std_b_insula = std(1000*b_insula);
b_Z_insula=mean_b_insula./std_b_insula;
b_P_insula = 2*normcdf(-1*abs(b_Z_insula),0,1);
%% NAc surface
figure;

pos_affect_pattern = tv;
pos_affect_pattern = apply_mask(pos_affect_pattern,nac);
st=tpaps([abs(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,1)-41),(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,2))]',pos_affect_pattern.dat');
    

points = fnplt(st);

surfqc(points{1}*2,points{2}*+pos_affect_pattern.volInfo.mat(2,2)+pos_affect_pattern.volInfo.mat(2,4),points{3})

 colormap(canlabCmap)
caxis([-.00008 .00003])
xlabel '|MNI_x|'
ylabel 'MNI_y'
zlabel 'Signature Coefficient'
drawnow; snapnow;

%% Insula surface
figure;
pos_affect_pattern = tv;
pos_affect_pattern = apply_mask(pos_affect_pattern,insula);
st=tpaps([abs(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,1)-41),(pos_affect_pattern.volInfo.xyzlist(~pos_affect_pattern.removed_voxels,2))]',pos_affect_pattern.dat');

points = fnplt(st);
surfqc(points{1}*2,points{2}*+pos_affect_pattern.volInfo.mat(2,2)+pos_affect_pattern.volInfo.mat(2,4),points{3})

colormap(canlabCmap)
xlabel '|MNI_x|'
ylabel 'MNI_y'
zlabel 'Signature Coefficient'

drawnow;snapnow;
basedir = fileparts(pwd);
datadir = [basedir filesep 'Data'];
resultsdir = [basedir filesep 'Results'];
codedir = [basedir filesep 'Code'];
anatdir = [basedir filesep 'AnatomicalMasks'];
tempdir = [codedir filesep 'intermediate_files'];

addpath(genpath(basedir));
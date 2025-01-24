% Create the directory that stores the simulations
directory = 'Sample_files';
mkdir(directory);

% End point of the simulation. Units in seconds.
tstop=4000;
% Number of realizations
random_seeds = 1;
% Number of Far1-GEF
n_FarGEF = 10;
% Minimum actin dissociation rate
mink = [0.005];

% Write config files
fid=fopen(sprintf('%s/run.sh',directory),'w');
for j = random_seeds
    for m = 1:numel(mink)
        curr_seed  = j;
        curr_fileprefix = sprintf('mink_%g-seeds_%d',mink(m),j);
        smoldyn_cfg(curr_fileprefix,directory,tstop,mink(m),curr_seed);
        cfg_name = sprintf('%s.cfg',curr_fileprefix);
    end
end
fclose(fid);
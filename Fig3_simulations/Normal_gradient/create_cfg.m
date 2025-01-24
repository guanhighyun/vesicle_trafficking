% Create the directory that stores the simulations
directory = 'Sample_files';
mkdir(directory);

% End point of the simulation. Units in seconds.
tstop=4000;
% Number of realizations
random_seeds = 1;
% Max actin dissociation rate in the center
maxk = [0.3];

% Create the script file for SLURM to run the jobs.
fid=fopen(sprintf('%s/run.sh',directory),'w');
for j = random_seeds
    for m = 1:numel(maxk)
        curr_seed  = j;
        curr_fileprefix = sprintf('maxk_%g-seeds_%d',maxk(m),j);
        smoldyn_cfg(curr_fileprefix,directory,tstop,maxk(m),curr_seed);
        cfg_name = sprintf('%s.cfg',curr_fileprefix);
    end
end
fclose(fid);
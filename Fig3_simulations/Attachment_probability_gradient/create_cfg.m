% Create the directory that stores the simulations
directory = 'Sample_files';
mkdir(directory);

% End point of the simulation. Units in seconds.
tstop=4000;
% Number of realizations
random_seeds = 1;
% Number of Far1-GEF
n_FarGEF = 30;
% Max actin attachment probability
maxP = [1];

fid=fopen(sprintf('%s/run.sh',directory),'w');
for j = random_seeds
    for m = 1:numel(maxP)
        curr_seed  = j;
        curr_fileprefix = sprintf('maxP_%g-seeds_%d',maxP(m),j);
        smoldyn_cfg(curr_fileprefix,directory,tstop,maxP(m),curr_seed);
        cfg_name = sprintf('%s.cfg',curr_fileprefix);
    end
end
fclose(fid);
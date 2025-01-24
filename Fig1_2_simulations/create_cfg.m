% Create the directory that stores the simulations
directory = 'Sample_files';
mkdir(directory);

% End point of the simulation. Units in seconds.
tstop=4000;
% Range of Bem1-GEF. Modify if needed.
n_BemGEF = 280; 
% Range of Cdc42. Modify if needed.
n_Cdc42 = 3000;
% Number of realizations
random_seeds = 1;
% Number of Far1-GEF
n_FarGEF = 30;
% Actin dissociation rate
k9b = 1/60;
% Number of actin cables
n_Actin = 10;

% Generate the configuration files
fid=fopen(sprintf('%s/run.sh',directory),'w');
for j = random_seeds
    for u = 1:numel(n_Cdc42)
        for v = 1:numel(n_BemGEF)
            for w = 1:numel(k9b)
                for z = 1:numel(n_Actin)
                    curr_seed  = j;
                    curr_fileprefix = sprintf('Cdc42_%g-Bem_%g-k9b_%g-n_actin_%d-seeds_%d',n_Cdc42(u),n_BemGEF(v),k9b(w),n_Actin(z),j);
                    smoldyn_cfg(curr_fileprefix,directory,tstop,n_Cdc42(u),n_BemGEF(v),k9b(w),n_Actin(z),curr_seed);
                    cfg_name = sprintf('%s.cfg',curr_fileprefix);
                end
            end
        end
    end
end
fclose(fid);
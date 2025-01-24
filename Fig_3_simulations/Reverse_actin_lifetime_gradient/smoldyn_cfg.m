% Adapted from "make_smoldyn_cfg.m" from Ramirez SA, Pablo M, Burk S,
% Lew DJ, Elston TC (2021). PLoS Comput Biol 17(7): e1008525.

% Inputs:
% fileprefix: the prefix of the cfg file.
% directory: directory where we write and store our cfg file. 
% tstop: time to stop the simulation
% samplingrate: Store the coordinates of molecules every 10 seconds. 
% Units in 0.1 ms.
% n_Cdc42: total number of Cdc42.
% n_BemGEF: total number of Bem1-GEF.
% mink: minimum Attached_actin dissociation rate.
% random_seed: a number which initializes a random number generator.

% Output:
% The cfg files will be generated in the "directory" folder.

function smoldyn_cfg(fileprefix,directory,tstop,mink,random_seed)

% configuration file for smoldyn
cfg_name = [directory '/' fileprefix '.cfg']; 
% xyz coords of species. they'll be written to the same directory as 
% the .cfg file.
xyz_name = [fileprefix '.xyz']; 
fid=fopen(cfg_name,'w');

fprintf(fid, 'random_seed %d\n',random_seed);
fprintf(fid, 'variable rho = 0.05\n');
fprintf(fid, 'variable rho_eps = 0.000010\n');
fprintf(fid, 'variable L = 8.8623\n');
fprintf(fid, 'variable k1a = 10\n');
fprintf(fid, 'variable k1b = 40\n');
fprintf(fid, 'variable P2a = 0.00053\n');
fprintf(fid, 'variable P2c = 0.018\n');
fprintf(fid, 'variable k2b = 0.35\n');
fprintf(fid, 'variable P3 = 0.018\n');
fprintf(fid, 'variable P4a = 0.00096\n');
fprintf(fid, 'variable k4b = 40\n');
fprintf(fid, 'variable k5a = 36\n');
fprintf(fid, 'variable k5b = 13\n');
fprintf(fid, 'variable P7  = 0.0256\n');
fprintf(fid, 'variable P6a  = 1\n');
fprintf(fid, 'variable P9a  = 0.04\n');
fprintf(fid, 'variable k9b  = %g\n',1/60);
fprintf(fid, 'variable P8a = 0.03\n');
fprintf(fid, 'variable k8b = 0.11\n');
fprintf(fid, 'variable P11 = 0.0013\n');
fprintf(fid, 'variable P12 = 2.5e-06\n');
fprintf(fid, 'variable k13 = 0.002\n');
fprintf(fid, 'variable k14 = 0.002\n');
fprintf(fid, 'variable k15 = 0.0004\n');

fprintf(fid, 'dim 2\n');
fprintf(fid, 'boundaries x -0.1 L+0.1\n');
fprintf(fid, 'boundaries y -0.1 L+0.1\n');
fprintf(fid, 'species Cdc42T Cdc42Dc Cdc42Dm BemGEFc BemGEFm BemGEF42 complex_Cdc42Dm_BemGEF42 complex_Cdc42Dm_BemGEFm\n');
fprintf(fid, 'species Vesicle Fused_vesicle complex_Cdc42T_Fused_vesicle Attached_actin Actin complex_Fused_vesicle_Attached_actin complex_Attached_actin_Cdc42T\n');

fprintf(fid, 'difc Cdc42T 0.0025\n');
fprintf(fid, 'difc Cdc42Dc 15\n');
fprintf(fid, 'difc Cdc42Dm 0.0025\n');
fprintf(fid, 'difc BemGEF42 0.0025\n');
fprintf(fid, 'difc BemGEFc 15\n');
fprintf(fid, 'difc BemGEFm 0.0025\n');
fprintf(fid, 'difc complex_Cdc42Dm_BemGEF42 0.0025\n');
fprintf(fid, 'difc complex_Cdc42Dm_BemGEFm 0.0025\n');
fprintf(fid, 'difc Vesicle 15\n');
fprintf(fid, 'difc Actin 15\n');
fprintf(fid, 'molecule_lists list1 list2 list3 list4 list5 list6 list7 list8 list9 list10 list11\n');
fprintf(fid, 'mol_list Cdc42T list1\n');
fprintf(fid, 'mol_list BemGEF42 list2\n');
fprintf(fid, 'mol_list Cdc42Dm list3\n');
fprintf(fid, 'mol_list Cdc42Dc list4\n');
fprintf(fid, 'mol_list BemGEFm list5\n');
fprintf(fid, 'mol_list BemGEFc list6\n');
fprintf(fid, 'mol_list Fused_vesicle list7\n');
fprintf(fid, 'mol_list Vesicle list8\n');
fprintf(fid, 'mol_list Attached_actin list9\n');
fprintf(fid, 'mol_list Actin list10\n');
fprintf(fid, 'mol_list A list11\n');
fprintf(fid, 'start_surface inner_walls\n');
fprintf(fid, 'action both all jump\n');
fprintf(fid, 'polygon both edge\n');
fprintf(fid, 'panel rect +x 0 0 L r1\n');
fprintf(fid, 'panel rect -x L 0 L r2\n');
fprintf(fid, 'panel rect +y 0 0 L r3\n');
fprintf(fid, 'panel rect -y 0 L L r4\n');
fprintf(fid, 'jump r1 front <-> r2 front\n');
fprintf(fid, 'jump r3 front <-> r4 front\n');
fprintf(fid, 'end_surface\n');

% Gaussian gradient of actin dissociation rate
L = 8.8623;
maxk = 1/60;
sigma = 0.5;
n_contour = 11;
sample_y = linspace(mink,maxk,n_contour);
sample_x = 1.4:-0.01:0;
f = (maxk-mink)*exp(-1/2*sample_x.^2/sigma^2)+mink;
k = nan(1,n_contour-1);
radius = nan(1,n_contour-1);
for i = 1:n_contour-1
    idx = (f>=sample_y(i)) & (f<sample_y(i+1));
    radius(i) = sample_x(find(f>=sample_y(i),1,'first'));
    k(i) = mean(f(idx));
end

radius = flip(radius);
for i = 1:n_contour-1
    fprintf(fid,'start_surface surf_%d\n',i);
    fprintf(fid,'polygon both edge\n');
    fprintf(fid,'action both all transmit\n');
    fprintf(fid,'panel sphere L/2 L/2 %g 100 panel_%d\n',radius(i),i);
    fprintf(fid,'end_surface\n');
end

for i = 1:n_contour-1
    fprintf(fid,'start_compartment cmpt_%d\n',i);
    if i ==1 
        fprintf(fid,'surface surf_%d\n',i);
        fprintf(fid,'point L/2 L/2\n');
    else
        fprintf(fid,'surface surf_%d\n',i);
        fprintf(fid,'surface surf_%d\n',i-1);
        fprintf(fid,'point L/2 %g\n',(L/2+radius(i-1)+L/2+radius(i))/2);
        fprintf(fid,'point L/2 %g\n',(L/2-radius(i-1)+L/2-radius(i))/2);
        fprintf(fid,'point %g L/2\n',(L/2+radius(i-1)+L/2+radius(i))/2);
        fprintf(fid,'point %g L/2\n',(L/2-radius(i-1)+L/2-radius(i))/2);
        fprintf(fid,'point %g %g\n',(L/2+radius(i-1)*sqrt(2)/2+L/2+radius(i)*sqrt(2)/2)/2,(L/2+radius(i-1)*sqrt(2)/2+L/2+radius(i)*sqrt(2)/2)/2);
        fprintf(fid,'point %g %g\n',(L/2+radius(i-1)*sqrt(2)/2+L/2+radius(i)*sqrt(2)/2)/2,(L/2-radius(i-1)*sqrt(2)/2+L/2-radius(i)*sqrt(2)/2)/2);
        fprintf(fid,'point %g %g\n',(L/2-radius(i-1)*sqrt(2)/2+L/2-radius(i)*sqrt(2)/2)/2,(L/2+radius(i-1)*sqrt(2)/2+L/2+radius(i)*sqrt(2)/2)/2);
        fprintf(fid,'point %g %g\n',(L/2-radius(i-1)*sqrt(2)/2+L/2-radius(i)*sqrt(2)/2)/2,(L/2-radius(i-1)*sqrt(2)/2+L/2-radius(i)*sqrt(2)/2)/2);
    end
    fprintf(fid,'end_compartment\n');
    fprintf(fid,'reaction compartment=cmpt_%d Attached_actin_%g Attached_actin -> Actin %g\n',i,i,k(i));
    fprintf(fid,'reaction surface=surf_%d Attached_actin_on_surf_%g Attached_actin -> Actin %g\n',i,i,k(i));
end

fprintf(fid,'start_compartment cmpt_%d\n',i+1);
fprintf(fid,'surface surf_%d\n',i);
fprintf(fid,'surface inner_walls\n');
fprintf(fid,'point L/2 L/2+%g\n',radius(end)+0.5);
fprintf(fid,'point L/2 L/2-%g\n',radius(end)+0.5);
fprintf(fid,'point L/2+%g L/2\n',radius(end)+0.5);
fprintf(fid,'point L/2-%g L/2\n',radius(end)+0.5);
fprintf(fid,'end_compartment\n');
fprintf(fid,'reaction compartment=cmpt_%d Attached_actin_%g Attached_actin -> Actin %g\n',i+1,i+1,maxk);
fprintf(fid,'reaction surface=inner_walls Attached_actin_on_inner_walls Attached_actin -> Actin %g\n',maxk);

for i = 1:n_contour-1
    if i ==1 
        fprintf(fid,'mol 1 A L/2 L/2\n');
    else
        fprintf(fid,'mol 1 A L/2 %g\n',(L/2+radius(i-1)+L/2+radius(i))/2);
        fprintf(fid,'mol 1 A L/2 %g\n',(L/2-radius(i-1)+L/2-radius(i))/2);
        fprintf(fid,'mol 1 A %g L/2\n',(L/2+radius(i-1)+L/2+radius(i))/2);
        fprintf(fid,'mol 1 A %g L/2\n',(L/2-radius(i-1)+L/2-radius(i))/2);
        fprintf(fid,'mol 1 A %g %g\n',(L/2+radius(i-1)*sqrt(2)/2+L/2+radius(i)*sqrt(2)/2)/2,(L/2+radius(i-1)*sqrt(2)/2+L/2+radius(i)*sqrt(2)/2)/2);
        fprintf(fid,'mol 1 A %g %g\n',(L/2+radius(i-1)*sqrt(2)/2+L/2+radius(i)*sqrt(2)/2)/2,(L/2-radius(i-1)*sqrt(2)/2+L/2-radius(i)*sqrt(2)/2)/2);
        fprintf(fid,'mol 1 A %g %g\n',(L/2-radius(i-1)*sqrt(2)/2+L/2-radius(i)*sqrt(2)/2)/2,(L/2+radius(i-1)*sqrt(2)/2+L/2+radius(i)*sqrt(2)/2)/2);
        fprintf(fid,'mol 1 A %g %g\n',(L/2-radius(i-1)*sqrt(2)/2+L/2-radius(i)*sqrt(2)/2)/2,(L/2-radius(i-1)*sqrt(2)/2+L/2-radius(i)*sqrt(2)/2)/2);
    end
end

fprintf(fid,'start_compartment full_domain\n');
fprintf(fid,'surface inner_walls\n');
fprintf(fid,'point L/2 L/2\n');
fprintf(fid,'end_compartment\n');

fprintf(fid, 'reaction BemGEF_cmtransition BemGEFc <-> BemGEFm k1a k1b\n');
fprintf(fid, 'reaction Cdc42_cmtransition Cdc42Dc <-> Cdc42Dm k5a k5b\n');
fprintf(fid, 'reaction Cdc42Dm_2_T_bindTo_BemGEFm Cdc42Dm + BemGEFm -> complex_Cdc42Dm_BemGEFm\n');
fprintf(fid, 'reaction Cdc42Dm_2_T_catBy_BemGEFm complex_Cdc42Dm_BemGEFm -> Cdc42T + BemGEFm\n');
fprintf(fid, 'reaction_probability Cdc42Dm_2_T_bindTo_BemGEFm P2a\n');
fprintf(fid, 'binding_radius Cdc42Dm_2_T_bindTo_BemGEFm rho\n');
fprintf(fid, 'reaction_probability Cdc42Dm_2_T_catBy_BemGEFm 1\n');
fprintf(fid, 'product_placement Cdc42Dm_2_T_catBy_BemGEFm unbindrad rho+rho_eps\n');
fprintf(fid, 'reaction Cdc42T_2_Cdc42Dm Cdc42T -> Cdc42Dm k2b\n');
fprintf(fid, 'reaction Cdc42Dm_2_T_bindTo_BemGEF42 Cdc42Dm + BemGEF42 -> complex_Cdc42Dm_BemGEF42\n');
fprintf(fid, 'reaction Cdc42Dm_2_T_catBy_BemGEF42 complex_Cdc42Dm_BemGEF42 -> Cdc42T + BemGEF42\n');
fprintf(fid, 'reaction_probability Cdc42Dm_2_T_bindTo_BemGEF42 P3\n');
fprintf(fid, 'binding_radius Cdc42Dm_2_T_bindTo_BemGEF42 rho\n');
fprintf(fid, 'reaction_probability Cdc42Dm_2_T_catBy_BemGEF42 1\n');
fprintf(fid, 'product_placement Cdc42Dm_2_T_catBy_BemGEF42 unbindrad rho+rho_eps\n');
fprintf(fid, 'reaction make_BemGEF42_fromm BemGEFm + Cdc42T <-> BemGEF42\n');
fprintf(fid, 'reaction_probability make_BemGEF42_frommfwd P4a\n');
fprintf(fid, 'binding_radius make_BemGEF42_frommfwd rho\n');
fprintf(fid, 'reaction_rate make_BemGEF42_frommrev k4b\n');
fprintf(fid, 'product_placement make_BemGEF42_frommrev unbindrad rho+rho_eps\n');
fprintf(fid, 'reaction make_BemGEF42_fromc BemGEFc + Cdc42T -> BemGEF42\n');
fprintf(fid, 'reaction_probability make_BemGEF42_fromc P7\n');
fprintf(fid, 'binding_radius make_BemGEF42_fromc rho\n');
fprintf(fid, 'reaction Attached_actin_polymerization Actin + Cdc42T -> complex_Attached_actin_Cdc42T\n');
fprintf(fid, 'reaction Attached_actin_polymerization_complex complex_Attached_actin_Cdc42T -> Attached_actin + Cdc42T\n');
fprintf(fid, 'reaction_probability Attached_actin_polymerization P9a\n');
fprintf(fid, 'reaction_probability Attached_actin_polymerization_complex 1\n');
fprintf(fid, 'binding_radius Attached_actin_polymerization rho\n');
fprintf(fid, 'product_placement Attached_actin_polymerization_complex unbindrad rho_eps\n');
fprintf(fid, 'reaction Yi_cyto2mem Vesicle + Attached_actin -> complex_Fused_vesicle_Attached_actin\n');
fprintf(fid, 'reaction Yi_cyto2mem_complex complex_Fused_vesicle_Attached_actin -> Fused_vesicle + Attached_actin\n');
fprintf(fid, 'reaction_probability Yi_cyto2mem P6a\n');
fprintf(fid, 'reaction_probability Yi_cyto2mem_complex 1\n');
fprintf(fid, 'binding_radius Yi_cyto2mem rho\n');
fprintf(fid, 'product_placement Yi_cyto2mem_complex unbindrad rho_eps\n');
fprintf(fid, 'reaction Yi_mem2cyto Fused_vesicle -> 0\n');
fprintf(fid, 'reaction_probability Yi_mem2cyto 1\n');
fprintf(fid, 'cmd I 2000000 80000000 1 longrangeforce Fused_vesicle Cdc42T 0 10^4 1e-05 3 (r^2+0.01)^0.5-r\n');
fprintf(fid, 'cmd I 2000000 80000000 1 longrangeforce Fused_vesicle BemGEF42 0 10^4 1e-05 3 (r^2+0.01)^0.5-r\n');
fprintf(fid, 'cmd I 2000000 80000000 1 longrangeforce Fused_vesicle BemGEFm 0 10^4 1e-05 3 (r^2+0.01)^0.5-r\n');
fprintf(fid, 'cmd I 2000000 80000000 1 longrangeforce Fused_vesicle Cdc42Dm 0 10^4 1e-05 3 (r^2+0.01)^0.5-r\n');
fprintf(fid, 'time_start 0\n');
fprintf(fid, 'time_stop %d\n',tstop);
fprintf(fid, 'time_step 0.0001\n');

fprintf(fid, 'compartment_mol 5 Vesicle full_domain\n');
fprintf(fid, 'cmd I 1 %d 60000 set compartment_mol 5 Vesicle full_domain\n',tstop/0.0001);

% Initial conditions
load('../initial_coordinates.mat')

for i = 1:height(Cdc42T)
    fprintf(fid, 'mol 1 Cdc42T %g %g\n',Cdc42T(i,1),Cdc42T(i,2));
end
for i = 1:height(BemGEF42)
    fprintf(fid, 'mol 1 BemGEF42 %g %g\n',BemGEF42(i,1),BemGEF42(i,2));
end
for i = 1:height(BemGEFc)
    fprintf(fid, 'mol 1 BemGEFc %g %g\n',BemGEFc(i,1),BemGEFc(i,2));
end
for i = 1:height(BemGEFm)
    fprintf(fid, 'mol 1 BemGEFm %g %g\n',BemGEFm(i,1),BemGEFm(i,2));
end
for i = 1:height(Cdc42Dc)
    fprintf(fid, 'mol 1 Cdc42Dc %g %g\n',Cdc42Dc(i,1),Cdc42Dc(i,2));
end
for i = 1:height(Cdc42Dm)
    fprintf(fid, 'mol 1 Cdc42Dm %g %g\n',Cdc42Dm(i,1),Cdc42Dm(i,2));
end
for i = 1:height(Attached_actin)
    fprintf(fid, 'mol 1 Attached_actin %g %g\n',Attached_actin(i,1),Attached_actin(i,2));
end


fprintf(fid, 'output_files %s\n', xyz_name);
fprintf(fid, 'cmd N 100000 molpos Cdc42T %s\n', xyz_name);
fprintf(fid, 'cmd N 100000 molpos BemGEF42 %s\n', xyz_name);
fprintf(fid, 'cmd N 100000 molpos Cdc42Dm %s\n', xyz_name);
fprintf(fid, 'cmd N 100000 molpos Cdc42Dc %s\n', xyz_name);
fprintf(fid, 'cmd N 100000 molpos BemGEFm %s\n', xyz_name);
fprintf(fid, 'cmd N 100000 molpos BemGEFc %s\n', xyz_name);
fprintf(fid, 'cmd N 100000 molpos Attached_actin %s\n', xyz_name);
fprintf(fid, 'cmd N 100000 molpos Actin %s\n', xyz_name);

end
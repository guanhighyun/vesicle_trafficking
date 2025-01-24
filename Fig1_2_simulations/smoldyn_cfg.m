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
% k9b: Attached_actin dissociation rate.
% n_Attached_actin: Number of Attached_actin cables
% random_seed: a number which initializes a random number generator.

% Output:
% The cfg files will be generated in the "directory" folder.

function smoldyn_cfg(fileprefix,directory,tstop,n_Cdc42,n_BemGEF,k9b,n_Attached_actin,random_seed)

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
fprintf(fid, 'variable P8a = 0.03\n');
fprintf(fid, 'variable k8b = 0.11\n');
fprintf(fid, 'variable P11 = 0.0013\n');
fprintf(fid, 'variable P12 = 2.5e-06\n');
fprintf(fid, 'variable k13 = 0.002\n');
fprintf(fid, 'variable k14 = 0.002\n');
fprintf(fid, 'variable k15 = 0.0004\n');

% Parameters for Attached_actin-mediated vesicle delivery
fprintf(fid, 'variable P6a  = 1\n');
fprintf(fid, 'variable k6b  = 10000\n');
fprintf(fid, 'variable P9a  = 0.04\n');
fprintf(fid, 'variable k9b  = %g\n',k9b);

fprintf(fid, 'dim 2\n');
fprintf(fid, 'boundaries x 0 L p\n');
fprintf(fid, 'boundaries y 0 L p\n');

fprintf(fid, 'species Cdc42T Cdc42Dc Cdc42Dm BemGEFc BemGEFm BemGEF42 complex_Cdc42Dm_BemGEF42 complex_Cdc42Dm_BemGEFm\n');
fprintf(fid, 'species Ram Rim Ric RaGEF FarGEF complex_Cdc42D_Ra complex_Ra_Cdc42T phe\n');
fprintf(fid, 'species Vesicle Fused_Vesicle complex_Cdc42T_Fused_Vesicle Attached_actin Actin complex_Fused_Vesicle_Attached_actin complex_Attached_actin_Cdc42T Attached_actin_Ram\n');
fprintf(fid, 'difc Cdc42T 0.0025\n');
fprintf(fid, 'difc Cdc42Dc 15\n');
fprintf(fid, 'difc Cdc42Dm 0.0025\n');
fprintf(fid, 'difc BemGEF42 0.0025\n');
fprintf(fid, 'difc BemGEFc 15\n');
fprintf(fid, 'difc BemGEFm 0.0025\n');
fprintf(fid, 'difc complex_Cdc42Dm_BemGEF42 0.0025\n');
fprintf(fid, 'difc complex_Cdc42Dm_BemGEFm 0.0025\n');
fprintf(fid, 'difc FarGEF 15\n');
fprintf(fid, 'difc Ram 0.0001\n');
fprintf(fid, 'difc RaGEF 0.0001\n');
fprintf(fid, 'difc Rim 0.0001\n');
fprintf(fid, 'difc Ric 15\n');
fprintf(fid, 'difc complex_Cdc42D_Ra 0.0001\n');
fprintf(fid, 'difc complex_Ra_Cdc42T 0.0001\n');
fprintf(fid, 'difc phe 150\n');
fprintf(fid, 'difc Vesicle 15\n');
fprintf(fid, 'difc Actin 15\n');

fprintf(fid, 'molecule_lists list1 list2 list3 list4 list5 list6 list7 list8 list9 list10 list11 list12 list13 list14 list15\n');
fprintf(fid, 'mol_list Cdc42T list1\n');
fprintf(fid, 'mol_list BemGEF42 list2\n');
fprintf(fid, 'mol_list Cdc42Dm list3\n');
fprintf(fid, 'mol_list Cdc42Dc list4\n');
fprintf(fid, 'mol_list BemGEFm list5\n');
fprintf(fid, 'mol_list BemGEFc list6\n');
fprintf(fid, 'mol_list Ram list7\n');
fprintf(fid, 'mol_list FarGEF list8\n');
fprintf(fid, 'mol_list Vesicle list9\n');
fprintf(fid, 'mol_list Attached_actin list10\n');
fprintf(fid, 'mol_list Actin list11\n');
fprintf(fid, 'mol_list Rim list12\n');
fprintf(fid, 'mol_list RaGEF list13\n');
fprintf(fid, 'mol_list Ric list14\n');
fprintf(fid, 'mol_list Fused_Vesicle list15\n');

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
fprintf(fid, 'reaction Rim_bindTo_phe  Rim + phe -> Ram\n');
fprintf(fid, 'reaction_probability Rim_bindTo_phe P11\n');
fprintf(fid, 'binding_radius Rim_bindTo_phe rho\n');
fprintf(fid, 'reaction Ric_recruitBy_Cdc42T Ric + Cdc42T -> complex_Ra_Cdc42T\n');
fprintf(fid, 'reaction Ric_2_Rim complex_Ra_Cdc42T -> Rim + Cdc42T\n');
fprintf(fid, 'reaction_probability Ric_recruitBy_Cdc42T P12\n');
fprintf(fid, 'reaction_probability Ric_2_Rim 1\n');
fprintf(fid, 'binding_radius Ric_recruitBy_Cdc42T rho\n');
fprintf(fid, 'product_placement Ric_2_Rim unbindrad rho+rho_eps\n');
fprintf(fid, 'reaction Ram_deactivation  Ram -> Rim + phe k13\n');
fprintf(fid, 'product_placement Ram_deactivation unbindrad rho+rho_eps\n');
fprintf(fid, 'reaction Ram_endocytosis Ram -> Ric + phe k14\n');
fprintf(fid, 'product_placement Ram_endocytosis unbindrad rho+rho_eps\n');
fprintf(fid, 'reaction Rim_endocytosis  Rim -> Ric k15\n');
fprintf(fid, 'reaction FarGEF_bindTo_Ra Ram + FarGEF -> RaGEF\n');
fprintf(fid, 'reaction_probability FarGEF_bindTo_Ra P8a\n');
fprintf(fid, 'binding_radius FarGEF_bindTo_Ra rho\n');
fprintf(fid, 'reaction RaGEF_dissociation RaGEF -> Ram + FarGEF k8b\n');
fprintf(fid, 'product_placement RaGEF_dissociation unbindrad rho+rho_eps\n');
fprintf(fid, 'reaction Cdc42D_2_T_bindTo_RaGEF  Cdc42Dm + RaGEF -> complex_Cdc42D_Ra\n');
fprintf(fid, 'reaction Cdc42D_2_T_catBy_RaGEF complex_Cdc42D_Ra -> Cdc42T + RaGEF\n');
fprintf(fid, 'reaction_probability Cdc42D_2_T_bindTo_RaGEF P2c\n');
fprintf(fid, 'binding_radius Cdc42D_2_T_bindTo_RaGEF rho\n');
fprintf(fid, 'reaction_probability Cdc42D_2_T_catBy_RaGEF 1\n');
fprintf(fid, 'product_placement Cdc42D_2_T_catBy_RaGEF unbindrad rho+rho_eps\n');
fprintf(fid, 'reaction phe_degradation  phe -> 0 6931.5\n');

% Reaction for Attached_actin-mediated vesicle delivery
fprintf(fid, 'reaction Attached_actin_polymerization Actin + Cdc42T -> complex_Attached_actin_Cdc42T\n');
fprintf(fid, 'reaction Attached_actin_polymerization_complex complex_Attached_actin_Cdc42T -> Attached_actin + Cdc42T\n');
fprintf(fid, 'reaction_probability Attached_actin_polymerization P9a\n');
fprintf(fid, 'reaction_probability Attached_actin_polymerization_complex 1\n');
fprintf(fid, 'binding_radius Attached_actin_polymerization rho\n');
fprintf(fid, 'product_placement Attached_actin_polymerization_complex unbindrad rho_eps\n');
fprintf(fid, 'reaction Attached_actin_depolymerization Attached_actin -> Actin k9b\n');
fprintf(fid, 'reaction V_cyto2mem Vesicle + Attached_actin -> complex_Fused_Vesicle_Attached_actin\n');
fprintf(fid, 'reaction V_cyto2mem_complex complex_Fused_Vesicle_Attached_actin -> Fused_Vesicle + Attached_actin\n');
fprintf(fid, 'reaction_probability V_cyto2mem P6a\n');
fprintf(fid, 'reaction_probability V_cyto2mem_complex 1\n');
fprintf(fid, 'binding_radius V_cyto2mem rho\n');
fprintf(fid, 'product_placement V_cyto2mem_complex unbindrad rho_eps\n');
fprintf(fid, 'reaction V_mem2cyto Fused_Vesicle -> 0\n');
fprintf(fid, 'reaction_probability V_mem2cyto 1\n');

% Vesicle fusion
fprintf(fid, 'cmd E longrangeforce Fused_Vesicle Cdc42T 0 10^4 1e-05 3 (r^2+0.01)^0.5-r\n');
fprintf(fid, 'cmd E longrangeforce Fused_Vesicle BemGEF42 0 10^4 1e-05 3 (r^2+0.01)^0.5-r\n');
fprintf(fid, 'cmd E longrangeforce Fused_Vesicle BemGEFm 0 10^4 1e-05 3 (r^2+0.01)^0.5-r\n');
fprintf(fid, 'cmd E longrangeforce Fused_Vesicle Cdc42Dm 0 10^4 1e-05 3 (r^2+0.01)^0.5-r\n');

fprintf(fid, 'time_start 0\n');
fprintf(fid, 'time_stop %g\n',tstop);
fprintf(fid, 'time_step 0.0001\n');
fprintf(fid, 'cmd @ 0.0001 gaussiansource Cdc42T %g L/2 0.2 L/2 0.2\n',n_Cdc42);
fprintf(fid, 'mol %g BemGEFc u u\n',n_BemGEF);

% Vesicle input
fprintf(fid, 'cmd @ 200 set mol 5 Vesicle u u\n');
fprintf(fid, 'cmd I 2000000 %d 60000 set mol 5 Vesicle u u\n',tstop/0.0001);

% Actin input
fprintf(fid, 'cmd @ 200 set mol %d Actin u u\n',n_Attached_actin);


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
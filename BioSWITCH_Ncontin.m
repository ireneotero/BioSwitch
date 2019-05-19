function out = BioSWITCH_Ncontin(input_file)
% BioSWITCH_Ncontin(input_file)  input_file 
% start a continuation analysis with Cl-Matcont for the reaction network and problem specifications defined in
% input_file
eval('run %s', input_file)




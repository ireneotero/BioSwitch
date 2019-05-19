function Results = BioSWITCH_Lpsearch(input_file)
% BioSWITCHS_Lpsearch(input_file)  input_file 
% search for a saddle node for the network and problem specifications defined in
% input_file


name_file=input_file;


eval(sprintf(name_file))

Results=MEIGO(problem,opts,'ESS');


eval(sprintf('save %s_Results Results',input_file))



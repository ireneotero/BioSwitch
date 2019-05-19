% Example 4

display_message('Evaluate Network');
BioSWITCH_Evaluate('EX4')

 display_message('Generate files for limit point search (system with mass conservation)');
 BioSWITCH_mkfiles_lpsearch_MassCon('EX4')
 
 display_message('Search for limit points (deficiency approach)');
 BioSWITCH_Lpsearch('EX4_input_lpsearch_MassCon_default')
 
 display_message('Generate files for continuation (system with mass conservation)');
 % indices to suppress from the odes
 BioSWITCH_mkfiles_ncontin_MassCon('EX4', [3 2 1]) % MUST BE IN ORDER
% 
% display_message('Check steady state and perform continuation (original network with conservation laws)');
 BioSWITCH_Ncontin('EX4_input_cont_MassCon_default')
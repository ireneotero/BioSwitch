% Example 3

display_message('Evaluate Network');
BioSWITCH_Evaluate('EX3')

display_message('Generate files for limit point search (system with mass conservation)');
BioSWITCH_mkfiles_lpsearch_MassCon('EX3')

display_message('Search for limit points (deficiency approach)');
BioSWITCH_Lpsearch('EX3_input_lpsearch_MassCon_default')

display_message('Generate files for continuation (system with mass conservation)');
% indices to suppress from the odes
BioSWITCH_mkfiles_ncontin_MassCon('EX3', [6 1 4]) % MUST BE IN ORDER

display_message('Check steady state and perform continuation (original network with conservation laws)');
BioSWITCH_Ncontin('EX3_input_cont_MassCon_default')
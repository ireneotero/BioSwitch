% Example G1S

display_message('Evaluate Network');
BioSWITCH_Evaluate('G1S')

display_message('Generate files for limit point search (system with mass conservation)');
display_message(' Please choose MANUAL option and vector number 63')

BioSWITCH_mkfiles_lpsearch_MassCon('G1S')

display_message('Search for limit points (deficiency approach)');
BioSWITCH_Lpsearch('G1S_input_lpsearch_MassCon_default')


display_message('Generate files for continuation (system with mass conservation)');
% indices to suppress from the odes
BioSWITCH_mkfiles_ncontin_MassCon('G1S', [7 3]) % MUST BE IN ORDER

display_message('Check steady state and perform continuation (original network with conservation laws)');
BioSWITCH_Ncontin('G1S_input_cont_MassCon_default')

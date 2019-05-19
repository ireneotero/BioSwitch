% Example 4b

display_message('Evaluate Network');
BioSWITCH_Evaluate('EX4')

display_message('Generate files for limit point search (semidiffusive system)');
BioSWITCH_mkfiles_lpsearch_SemiDiff('EX4', [1, 2, 3])

display_message('Search for limit points (injectivity approach)');
BioSWITCH_Lpsearch('EX4_input_lpsearch_SemiDiff_default')

display_message('Generate files for continuation (semidiffusive system)');
% indices to suppress from the odes
BioSWITCH_mkfiles_ncontin_SemiDiff('EX4', [3,2,1]) % MUST BE IN ORDER

display_message('Check Steady state and Perform continuation (semidiffusive network)');
BioSWITCH_Ncontin('EX4_input_cont_SemiDiff_default')
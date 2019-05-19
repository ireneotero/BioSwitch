% Example 3b

display_message('Evaluate Network');
BioSWITCH_Evaluate('EX3')

display_message('Generate files for limit point search (semidiffusive system)');
BioSWITCH_mkfiles_lpsearch_SemiDiff('EX3', [1, 2, 6])

display_message('Search for limit points (injectivity approach)');
BioSWITCH_Lpsearch('EX3_input_lpsearch_SemiDiff_default')

display_message('Generate files for continuation (semidiffusive system)');
% indices to suppress from the odes
BioSWITCH_mkfiles_ncontin_SemiDiff('EX3',[6,1,2]) % MUST BE IN ORDER

display_message('Check Steady state and Perform continuation (semidiffusive network)');
BioSWITCH_Ncontin('EX3_input_cont_SemiDiff_default')
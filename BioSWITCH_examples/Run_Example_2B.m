% Example 2b

display_message('Evaluate Network');
BioSWITCH_Evaluate('EX2tr')

display_message('Generate files for limit point search (semidiffusive system)');
BioSWITCH_mkfiles_lpsearch_SemiDiff('EX2tr', [1, 5])

display_message('Search for limit points (injectivity approach)');
BioSWITCH_Lpsearch('EX2tr_input_lpsearch_SemiDiff_default')


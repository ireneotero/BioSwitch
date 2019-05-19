% Example 0

display_message('Evaluate Network');
BioSWITCH_Evaluate('EX0')

display_message('Generate files for limit point search (semidiffusive system)');
BioSWITCH_mkfiles_lpsearch_SemiDiff('EX0', [1 2])

display_message('Search for limit points (injectivity approach)');
Results = BioSWITCH_Lpsearch('EX0_input_lpsearch_SemiDiff_default')

if Results.fbest > 1e-7    
    disp('No Limit Pount found, we start now an equilibrium continuation from the optimum')
    pause(5)
end

display_message('Generate files for continuation (semidiffusive system)');
BioSWITCH_mkfiles_ncontin_SemiDiff('EX0', [1 2])

display_message('Check Steady state and Perform continuation');
BioSWITCH_Ncontin('EX0_input_cont_SemiDiff_default')

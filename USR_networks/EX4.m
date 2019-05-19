function [ESP, COM, KNS] = EX4
% Reference: model C3dd from J. Saez-Rodriguez et al, 2005,  IEEE Proceedings Systems Biology

ESP={'A','E1','E2','Ap','AE1','ApE2','AAp'};
COM={'A+E1','AE1','Ap+E1','Ap+E2','ApE2','A+E2','A+Ap','AAp','2*Ap'}; 
KNS={'k_1_2' 'k_2_1' 'k_2_3' 'k_4_5' 'k_5_4' 'k_5_6' 'k_7_8' 'k_8_7' 'k_8_9'};

end
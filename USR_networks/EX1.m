function [ESP, COM, KNS] = EX1
% Toggle switch; bistable
% Reference: Otero-Muras, Yordanov, Stelling BMC Syst Biol 2014

ESP={'G1','P1','G2','G1P2','G2P1','G2P1P1','P2'};
COM={'G1','G1+P1','G2','G2+P2','G1+P2','G1P2','G2+P1','G2P1','G2P1+P1','G2P1P1','P1','P2','0' }; 
KNS={'k_1_2' 'k_3_4' 'k_5_6' 'k_6_5' 'k_7_8' 'k_8_7'  'k_9_10' 'k_10_9' 'k_11_13' 'k_12_13'  };

end
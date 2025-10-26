function params = init_parameters()

    params.fc = 0.1e12;
    params.BW_narrow = 100e6;
    params.BW_wide = 10e9;
    

    params.N = 128;
    params.NRF = 4;
    params.K = 4;
    params.M = 6*6;
    params.NT = 8; % TTD per RF chain
    params.Mc = 10; % subcarriers
    params.LCP = 4;
    
 
    params.Pt_dBm = -20:5:50;
    params.Pt = 10.^(params.Pt_dBm/10)/1000;
    
  
    params.PBS = 3;
    params.PBB = 0.3;
    params.PRF = 0.2;
    params.PPS = 0.03;
    params.PTTD = 0.1;
    params.PUE = 0.1;
    params.xi = 0.1;
    
    params.NumIter = 100;
end

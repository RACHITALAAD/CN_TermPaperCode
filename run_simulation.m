function [SE_HB, SE_TTD, EE_HB, EE_TTD] = run_simulation(params)

L = 4;      % BS-STARS paths
Lk = 4;     % STARS-User paths

SE_HB = zeros(length(params.Pt),1);
SE_TTD = zeros(length(params.Pt),1);
EE_HB = zeros(length(params.Pt),1);
EE_TTD = zeros(length(params.Pt),1);

for p = 1:length(params.Pt)
    
    se_temp_HB = zeros(params.NumIter,1);
    se_temp_TTD = zeros(params.NumIter,1);
    
    for iter = 1:params.NumIter
        % BS to STARS
        H_BS_STARS = sum((randn(params.N,params.M,L)+1j*randn(params.N,params.M,L))/sqrt(2),3);
        
        % STARS to Users
        H_STARS_UE = zeros(params.M, params.K);
        for k_idx = 1:params.K
            H_STARS_UE(:,k_idx) = sum((randn(params.M,1,Lk)+1j*randn(params.M,1,Lk))/sqrt(2),3);
        end
        
 
        theta = exp(1j*2*pi*rand(params.M,1)); % Independent phase shifts
        
    
        F_PS = exp(1j*2*pi*rand(params.N,params.NRF)); % Analog
        F_BB = sqrt(params.Pt(p)/params.NRF)*(randn(params.NRF,params.K)+1j*randn(params.NRF,params.K))/sqrt(2);
        F_eff = F_PS*F_BB; % Effective hybrid beamformer
        
        % Compute SE for narrowband HB
        Y_HB = zeros(params.K,1);
        for k_idx = 1:params.K
            % Element-wise multiplication + matrix multiplication fix
            h_eff = H_STARS_UE(:,k_idx).' .* theta.'; % 1xM
            h_eff = h_eff * H_BS_STARS.';             % 1xN
            signal = h_eff * F_eff(:,k_idx);
            interference = sum(h_eff*F_eff,2) - signal;
            noise = sqrt(1e-9)*(randn+1j*randn);
            SINR = abs(signal)^2/(abs(interference)^2 + abs(noise)^2);
            Y_HB(k_idx) = log2(1+SINR);
        end
        se_temp_HB(iter) = sum(Y_HB);
     
        SE_sub = zeros(params.Mc,1);
        for m = 1:params.Mc
            f = params.fc + (m-1)*params.BW_wide/params.Mc;
            
            H_BS_STARS_f = H_BS_STARS .* exp(-1j*2*pi*(f-params.fc)*rand(params.N,params.M));
            H_STARS_UE_f = H_STARS_UE .* exp(-1j*2*pi*(f-params.fc)*rand(params.M,params.K));
            
            F_TTD = exp(1j*2*pi*rand(params.N,params.NRF)); % Analog
            F_BB_f = sqrt(params.Pt(p)/params.NRF)*(randn(params.NRF,params.K)+1j*randn(params.NRF,params.K))/sqrt(2);
            F_eff_f = F_TTD*F_BB_f;
            
            % Compute SE per subcarrier
            Y_sub = zeros(params.K,1);
            for k_idx = 1:params.K
                h_eff_f = H_STARS_UE_f(:,k_idx).' .* theta.'; % 1xM
                h_eff_f = h_eff_f * H_BS_STARS_f.';            % 1xN
                signal = h_eff_f * F_eff_f(:,k_idx);
                interference = sum(h_eff_f*F_eff_f,2) - signal;
                noise = sqrt(1e-9)*(randn+1j*randn);
                SINR = abs(signal)^2/(abs(interference)^2 + abs(noise)^2);
                Y_sub(k_idx) = log2(1+SINR);
            end
            SE_sub(m) = sum(Y_sub);
        end
        se_temp_TTD(iter) = mean(SE_sub); % Average over subcarriers
        
    end
    
   
    SE_HB(p) = mean(se_temp_HB);
    SE_TTD(p) = mean(se_temp_TTD);
    
    P_HB = params.PBS + params.PBB + params.NRF*params.PRF + params.NRF*params.PPS + params.xi*SE_HB(p);
    P_TTD = params.PBS + params.PBB + params.NRF*params.PRF + params.NRF*params.NT*params.PTTD + params.NRF*params.PPS + params.xi*SE_TTD(p);
    
    EE_HB(p) = SE_HB(p)/P_HB;
    EE_TTD(p) = SE_TTD(p)/P_TTD;
    
end
end

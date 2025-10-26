function performance_vs_STARS_elements(params)

M_list = [4, 6, 8, 10]; % Grid size (MxM elements)
SE_HB_M = zeros(length(M_list),1);
SE_TTD_M = zeros(length(M_list),1);
EE_HB_M = zeros(length(M_list),1);
EE_TTD_M = zeros(length(M_list),1);

SE_HB_noSTARS_M = zeros(length(M_list),1);
SE_TTD_noSTARS_M = zeros(length(M_list),1);
EE_HB_noSTARS_M = zeros(length(M_list),1);
EE_TTD_noSTARS_M = zeros(length(M_list),1);

% Find closest transmit power to 20 dBm
[~, Pt_idx] = min(abs(params.Pt_dBm - 20));
Pt = params.Pt(Pt_idx);

for m_idx = 1:length(M_list)
    M = M_list(m_idx)^2; % Total elements
    se_temp_HB = zeros(params.NumIter,1);
    se_temp_TTD = zeros(params.NumIter,1);
    se_temp_HB_noSTARS = zeros(params.NumIter,1);
    se_temp_TTD_noSTARS = zeros(params.NumIter,1);
    
    for iter = 1:params.NumIter
        L = 4; Lk = 4;

        % With STARS
        H_BS_STARS = sum((randn(params.N,M,L)+1j*randn(params.N,M,L))/sqrt(2),3);  % N x M
        H_STARS_UE = zeros(M, params.K);
        for k_idx = 1:params.K
            H_STARS_UE(:,k_idx) = sum((randn(M,1,Lk)+1j*randn(M,1,Lk))/sqrt(2),3);
        end
        theta = exp(1j*2*pi*rand(M,1)); % M x 1

        F_PS = exp(1j*2*pi*rand(params.N,params.NRF)); % N x NRF
        F_BB = sqrt(Pt/params.NRF)*(randn(params.NRF,params.K)+1j*randn(params.NRF,params.K))/sqrt(2); % NRF x K
        F_eff = F_PS*F_BB; % N x K

        Y_HB = zeros(params.K,1);
        Y_HB_noSTARS = zeros(params.K,1);
        for k_idx = 1:params.K
            % With STARS
            h_eff = (theta.' .* H_STARS_UE(:,k_idx).') * H_BS_STARS.'; % 1 x N
            signal = h_eff * F_eff(:,k_idx);
            interference = sum(h_eff*F_eff,2) - signal;
            noise = sqrt(1e-9)*(randn+1j*randn);
            SINR = abs(signal)^2/(abs(interference)^2 + abs(noise)^2);
            Y_HB(k_idx) = log2(1+SINR);

            % Without STARS (direct BS to User)
            h_eff_noSTARS = randn(1, params.N) + 1j*randn(1, params.N);
            signal_noSTARS = h_eff_noSTARS * F_eff(:,k_idx);
            interference_noSTARS = sum(h_eff_noSTARS*F_eff,2) - signal_noSTARS;
            SINR_noSTARS = abs(signal_noSTARS)^2 / (abs(interference_noSTARS)^2 + abs(noise)^2);
            Y_HB_noSTARS(k_idx) = log2(1+SINR_noSTARS);
        end
        se_temp_HB(iter) = sum(Y_HB);
        se_temp_HB_noSTARS(iter) = sum(Y_HB_noSTARS);

        % Wideband TTD 
        SE_sub = zeros(params.Mc,1);
        SE_sub_noSTARS = zeros(params.Mc,1);
        for sub = 1:params.Mc
            f = params.fc + (sub-1)*params.BW_wide/params.Mc;
            H_BS_STARS_f = H_BS_STARS .* exp(-1j*2*pi*(f-params.fc)*rand(params.N,M));
            H_STARS_UE_f = H_STARS_UE .* exp(-1j*2*pi*(f-params.fc)*rand(M,params.K));

            F_TTD = exp(1j*2*pi*rand(params.N,params.NRF));
            F_BB_f = sqrt(Pt/params.NRF)*(randn(params.NRF,params.K)+1j*randn(params.NRF,params.K))/sqrt(2);
            F_eff_f = F_TTD*F_BB_f;

            Y_sub = zeros(params.K,1);
            Y_sub_noSTARS = zeros(params.K,1);
            for k_idx = 1:params.K
                % With STARS
                h_eff_f = (theta.' .* H_STARS_UE_f(:,k_idx).') * H_BS_STARS_f.'; % 1 x N
                signal = h_eff_f * F_eff_f(:,k_idx);
                interference = sum(h_eff_f*F_eff_f,2) - signal;
                noise = sqrt(1e-9)*(randn+1j*randn);
                SINR = abs(signal)^2/(abs(interference)^2 + abs(noise)^2);
                Y_sub(k_idx) = log2(1+SINR);

                % Without STARS
                h_eff_noSTARS_f = randn(1, params.N) + 1j*randn(1, params.N);
                signal_noSTARS = h_eff_noSTARS_f * F_eff_f(:,k_idx);
                interference_noSTARS = sum(h_eff_noSTARS_f*F_eff_f,2) - signal_noSTARS;
                SINR_noSTARS = abs(signal_noSTARS)^2 / (abs(interference_noSTARS)^2 + abs(noise)^2);
                Y_sub_noSTARS(k_idx) = log2(1+SINR_noSTARS);
            end
            SE_sub(sub) = sum(Y_sub);
            SE_sub_noSTARS(sub) = sum(Y_sub_noSTARS);
        end
        se_temp_TTD(iter) = mean(SE_sub);
        se_temp_TTD_noSTARS(iter) = mean(SE_sub_noSTARS);
    end

    % Average Monte Carlo
    SE_HB_M(m_idx) = mean(se_temp_HB);
    SE_TTD_M(m_idx) = mean(se_temp_TTD);
    SE_HB_noSTARS_M(m_idx) = mean(se_temp_HB_noSTARS);
    SE_TTD_noSTARS_M(m_idx) = mean(se_temp_TTD_noSTARS);

    % Energy Efficiency
    P_HB = params.PBS + params.PBB + params.NRF*params.PRF + params.NRF*params.PPS + params.xi*SE_HB_M(m_idx);
    P_TTD = params.PBS + params.PBB + params.NRF*params.PRF + params.NRF*params.NT*params.PTTD + params.NRF*params.PPS + params.xi*SE_TTD_M(m_idx);
    P_HB_noSTARS = params.PBS + params.PBB + params.NRF*params.PRF + params.NRF*params.PPS + params.xi*SE_HB_noSTARS_M(m_idx);
    P_TTD_noSTARS = params.PBS + params.PBB + params.NRF*params.PRF + params.NRF*params.NT*params.PTTD + params.NRF*params.PPS + params.xi*SE_TTD_noSTARS_M(m_idx);

    EE_HB_M(m_idx) = SE_HB_M(m_idx)/P_HB;
    EE_TTD_M(m_idx) = SE_TTD_M(m_idx)/P_TTD;
    EE_HB_noSTARS_M(m_idx) = SE_HB_noSTARS_M(m_idx)/P_HB_noSTARS;
    EE_TTD_noSTARS_M(m_idx) = SE_TTD_noSTARS_M(m_idx)/P_TTD_noSTARS;

    % Print results 
    fprintf('M = %d elements:\n', M);
    fprintf('  Without STARS: SE_HB = %.4f, SE_TTD = %.4f, EE_HB = %.4f, EE_TTD = %.4f\n', ...
        SE_HB_noSTARS_M(m_idx), SE_TTD_noSTARS_M(m_idx), EE_HB_noSTARS_M(m_idx), EE_TTD_noSTARS_M(m_idx));
    fprintf('  With STARS   : SE_HB = %.4f, SE_TTD = %.4f, EE_HB = %.4f, EE_TTD = %.4f\n\n', ...
        SE_HB_M(m_idx), SE_TTD_M(m_idx), EE_HB_M(m_idx), EE_TTD_M(m_idx));
end

% Plot Spectral Efficiency
figure;
plot(M_list.^2, SE_HB_noSTARS_M,'--o','LineWidth',2); hold on;
plot(M_list.^2, SE_TTD_noSTARS_M,'--s','LineWidth',2);
plot(M_list.^2, SE_HB_M,'-o','LineWidth',2);
plot(M_list.^2, SE_TTD_M,'-s','LineWidth',2);
xlabel('Number of STARS Elements M');
ylabel('Spectral Efficiency (bits/s/Hz)');
legend('HB No STARS','TTD No STARS','HB With STARS','TTD With STARS');
grid on;
title('Spectral Efficiency vs Number of STARS Elements');

% Plot Energy Efficiency
figure;
plot(M_list.^2, EE_HB_noSTARS_M,'--o','LineWidth',2); hold on;
plot(M_list.^2, EE_TTD_noSTARS_M,'--s','LineWidth',2);
plot(M_list.^2, EE_HB_M,'-o','LineWidth',2);
plot(M_list.^2, EE_TTD_M,'-s','LineWidth',2);
xlabel('Number of STARS Elements M');
ylabel('Energy Efficiency (bits/Joule)');
legend('HB No STARS','TTD No STARS','HB With STARS','TTD With STARS');
grid on;
title('Energy Efficiency vs Number of STARS Elements');

end

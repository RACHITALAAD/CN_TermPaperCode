function plot_results(params, SE_HB, SE_TTD, EE_HB, EE_TTD)
    figure;
    plot(params.Pt_dBm, SE_HB,'-o','LineWidth',2); hold on;
    plot(params.Pt_dBm, SE_TTD,'-s','LineWidth',2);
    xlabel('Transmit Power Pt (dBm)');
    ylabel('Spectral Efficiency (bits/s/Hz)');
    legend('Hybrid Beamforming','TTD-based Hybrid');
    grid on; title('Spectral Efficiency vs Transmit Power');
    
    figure;
    plot(params.Pt_dBm, EE_HB,'-o','LineWidth',2); hold on;
    plot(params.Pt_dBm, EE_TTD,'-s','LineWidth',2);
    xlabel('Transmit Power Pt (dBm)');
    ylabel('Energy Efficiency (bits/Joule)');
    legend('Hybrid Beamforming','TTD-based Hybrid');
    grid on; title('Energy Efficiency vs Transmit Power');
end

clc; clear; close all;

params = init_parameters();

[SE_HB, SE_TTD, EE_HB, EE_TTD] = run_simulation(params);

plot_results(params, SE_HB, SE_TTD, EE_HB, EE_TTD);

performance_vs_STARS_elements(params);

disp('All simulations completed!');

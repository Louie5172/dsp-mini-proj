% this is the script for the higher level stuff:

time_domain_sim;

energy = get_energy(spatial_sin2, fs);

f_sin_mix = 80;
[disp1, press1] = get_sin_response(f_sin_mix, fs, coeff, n, Gp_idx, b1_idx, b2_idx, Gz_idx, a1_idx, a2_idx);
[disp2, press2] = get_sin_response(f_sin, fs, coeff, n, Gp_idx, b1_idx, b2_idx, Gz_idx, a1_idx, a2_idx);

displacement_sin_mix = disp1 + disp2;
pressure_sin_mix = press1 + press2;

plot_at_time(displacement_sin_mix, 1000);

fprintf("debug??");


function v_matrix = get_energy(displacement, fs)
    fc = 30000; % rate at which the switch is closed
    c_0 = exp(-2 * pi * fc / fs);

    e_matrix = max(displacement, 0);
    v_matrix = (1 - c_0) * e_matrix;

    % we go across and apply the same func for each row:
    for i = 2:size(displacement, 1)
        v_matrix(i, :) = v_matrix(i, :) + c_0 * v_matrix(i - 1, :);
    end
end

% plot displacement across filters at a given time
function plot_at_time(displacement, sample)
    figure;
    plot(1:size(displacement, 2), displacement(sample, :));

    xlabel('Filter Number');
    ylabel('Displacement');
    title(sprintf('Displacement at Sample %d', sample));
end
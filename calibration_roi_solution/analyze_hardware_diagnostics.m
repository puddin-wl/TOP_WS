function report = analyze_hardware_diagnostics(run_dir)
% analyze_hardware_diagnostics Quantify one hardware_diagnostics run folder.
%
% Usage:
%   analyze_hardware_diagnostics
%   analyze_hardware_diagnostics('E:\...\hardware_diagnostics\run_YYYYMMDD_HHMMSS')

    if nargin < 1 || isempty(run_dir)
        base_dir = fullfile(fileparts(mfilename('fullpath')), 'hardware_diagnostics');
        runs = dir(fullfile(base_dir, 'run_*'));
        if isempty(runs)
            error('No hardware diagnostic runs found under %s.', base_dir);
        end
        [~, newest_idx] = max([runs.datenum]);
        run_dir = fullfile(runs(newest_idx).folder, runs(newest_idx).name);
    end

    diag_file = fullfile(run_dir, 'hardware_diagnostics.mat');
    if ~exist(diag_file, 'file')
        error('Missing diagnostics file: %s', diag_file);
    end

    S = load(diag_file);
    diag_data = S.diag;
    loops = 1:diag_data.params.closed_loop_iters;

    target_amp = diag_data.Target_Amp_ROI;
    target_level = mean(target_amp(:));
    target_is_flat = max(abs(target_amp(:) - target_level)) < 1e-9;

    roi_raw = diag_data.sampled_roi_smooth_stack;
    amp_norm = diag_data.normalized_amp_roi_stack;
    amp_err = amp_norm - target_amp;
    has_point_sampling = isfield(diag_data, 'point_sampled_roi_smooth_stack') && ...
        isfield(diag_data, 'point_normalized_amp_roi_stack');
    if has_point_sampling
        point_roi_raw = diag_data.point_sampled_roi_smooth_stack;
        point_amp_err = diag_data.point_normalized_amp_roi_stack - target_amp;
    end

    report.run_dir = run_dir;
    report.feedback_sampling = 'area_average';
    if isfield(diag_data.params, 'roi_area_sample_count')
        report.roi_area_sample_count = diag_data.params.roi_area_sample_count;
    else
        report.feedback_sampling = 'point_or_legacy';
        report.roi_area_sample_count = 1;
    end
    report.best_loop = diag_data.best_loop;
    report.best_rms = diag_data.best_err_rms;
    report.first_rms = diag_data.err_rms(1);
    report.last_rms = diag_data.err_rms(end);
    report.rms_drop_pct = 100 * (report.first_rms - report.last_rms) / report.first_rms;
    report.target_amp_mean = target_level;
    report.target_is_flat = target_is_flat;
    report.sample_oob_max = max(diag_data.sample_oob_ratio);
    report.sample_valid_min = min(diag_data.sample_valid_ratio);

    roi_mean = squeeze(mean(roi_raw, [1 2]));
    roi_std = squeeze(std(reshape(roi_raw, [], size(roi_raw, 3)), 0, 1));
    roi_min = squeeze(min(reshape(roi_raw, [], size(roi_raw, 3)), [], 1));
    roi_max = squeeze(max(reshape(roi_raw, [], size(roi_raw, 3)), [], 1));
    report.intensity_mean_first_last = [roi_mean(1), roi_mean(end)];
    report.intensity_std_first_last = [roi_std(1), roi_std(end)];
    report.intensity_cv_first_last = [roi_std(1) / roi_mean(1), roi_std(end) / roi_mean(end)];
    report.intensity_min_first_last = [roi_min(1), roi_min(end)];
    report.intensity_max_first_last = [roi_max(1), roi_max(end)];

    best_err = amp_err(:, :, report.best_loop);
    best_raw = roi_raw(:, :, report.best_loop);
    report.best_error_mean = mean(best_err(:));
    report.best_error_std = std(best_err(:));
    report.best_error_min = min(best_err(:));
    report.best_error_max = max(best_err(:));
    report.best_error_abs_p95 = prctile(abs(best_err(:)), 95);
    report.best_intensity_mean = mean(best_raw(:));
    report.best_intensity_std = std(best_raw(:));
    report.best_intensity_cv = report.best_intensity_std / report.best_intensity_mean;
    report.best_intensity_min = min(best_raw(:));
    report.best_intensity_max = max(best_raw(:));

    if has_point_sampling
        point_best_err = point_amp_err(:, :, report.best_loop);
        point_best_raw = point_roi_raw(:, :, report.best_loop);
        report.point_best_error_rms = sqrt(mean(point_best_err(:).^2));
        report.point_best_error_abs_p95 = prctile(abs(point_best_err(:)), 95);
        report.point_best_intensity_cv = std(point_best_raw(:)) / mean(point_best_raw(:));
        report.point_area_rms_ratio = report.point_best_error_rms / report.best_rms;
    else
        point_best_err = [];
        report.point_best_error_rms = NaN;
        report.point_best_error_abs_p95 = NaN;
        report.point_best_intensity_cv = NaN;
        report.point_area_rms_ratio = NaN;
    end

    h = size(best_err, 1);
    w = size(best_err, 2);
    edge_mask = false(h, w);
    edge_width = max(1, round(min(h, w) * 0.15));
    edge_mask(1:edge_width, :) = true;
    edge_mask(end-edge_width+1:end, :) = true;
    edge_mask(:, 1:edge_width) = true;
    edge_mask(:, end-edge_width+1:end) = true;
    center_mask = ~edge_mask;
    report.edge_abs_error_mean = mean(abs(best_err(edge_mask)));
    report.center_abs_error_mean = mean(abs(best_err(center_mask)));
    report.edge_center_abs_error_ratio = report.edge_abs_error_mean / report.center_abs_error_mean;

    [xg, yg] = meshgrid(1:w, 1:h);
    X = [ones(numel(xg), 1), xg(:), yg(:)];
    plane_coef = X \ best_err(:);
    plane_fit = reshape(X * plane_coef, h, w);
    report.error_plane_coeff = plane_coef(:).';
    report.error_plane_rms = sqrt(mean(plane_fit(:).^2));
    report.error_residual_rms = sqrt(mean((best_err(:) - plane_fit(:)).^2));
    report.error_plane_energy_fraction = report.error_plane_rms^2 / max(mean(best_err(:).^2), eps);

    err_flat = abs(best_err(:));
    [~, order] = sort(err_flat, 'descend');
    top_n = min(10, numel(order));
    [top_r, top_c] = ind2sub([h, w], order(1:top_n));
    report.top_abs_error_pixels = table(top_r, top_c, best_err(order(1:top_n)), ...
        best_raw(order(1:top_n)), ...
        'VariableNames', {'row', 'col', 'amp_error', 'intensity'});

    row_mean = mean(best_err, 2);
    col_mean = mean(best_err, 1);
    report.best_error_row_mean = row_mean;
    report.best_error_col_mean = col_mean;

    out_txt = fullfile(run_dir, 'analysis_report.txt');
    fid = fopen(out_txt, 'w');
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, 'Hardware diagnostics analysis\n');
    fprintf(fid, 'Run: %s\n\n', run_dir);
    fprintf(fid, 'RMS: first %.6f, last %.6f, best loop %d = %.6f, drop %.2f%%\n', ...
        report.first_rms, report.last_rms, report.best_loop, report.best_rms, report.rms_drop_pct);
    fprintf(fid, 'Feedback sampling: %s, area sample count = %d x %d\n', ...
        report.feedback_sampling, report.roi_area_sample_count, report.roi_area_sample_count);
    if has_point_sampling
        fprintf(fid, 'Point-sampling comparison at best loop: RMS %.6f, abs p95 %.6f, CV %.6f, point/area RMS ratio %.6f\n', ...
            report.point_best_error_rms, report.point_best_error_abs_p95, ...
            report.point_best_intensity_cv, report.point_area_rms_ratio);
    end
    fprintf(fid, 'Sampling: max OOB %.8f, min valid %.8f\n', ...
        report.sample_oob_max, report.sample_valid_min);
    fprintf(fid, 'Target amplitude: mean %.6f, flat=%d\n\n', ...
        report.target_amp_mean, report.target_is_flat);
    fprintf(fid, 'Intensity ROI first/last:\n');
    fprintf(fid, '  mean %.6f -> %.6f\n', report.intensity_mean_first_last);
    fprintf(fid, '  std  %.6f -> %.6f\n', report.intensity_std_first_last);
    fprintf(fid, '  CV   %.6f -> %.6f\n', report.intensity_cv_first_last);
    fprintf(fid, '  min  %.6f -> %.6f\n', report.intensity_min_first_last);
    fprintf(fid, '  max  %.6f -> %.6f\n\n', report.intensity_max_first_last);
    fprintf(fid, 'Best loop spatial error:\n');
    fprintf(fid, '  amp error mean/std/min/max = %.6f / %.6f / %.6f / %.6f\n', ...
        report.best_error_mean, report.best_error_std, report.best_error_min, report.best_error_max);
    fprintf(fid, '  abs error p95 = %.6f\n', report.best_error_abs_p95);
    fprintf(fid, '  intensity mean/std/CV/min/max = %.6f / %.6f / %.6f / %.6f / %.6f\n', ...
        report.best_intensity_mean, report.best_intensity_std, report.best_intensity_cv, ...
        report.best_intensity_min, report.best_intensity_max);
    fprintf(fid, '  edge abs error %.6f, center abs error %.6f, ratio %.6f\n', ...
        report.edge_abs_error_mean, report.center_abs_error_mean, report.edge_center_abs_error_ratio);
    fprintf(fid, '  plane energy fraction %.6f, plane RMS %.6f, residual RMS %.6f\n\n', ...
        report.error_plane_energy_fraction, report.error_plane_rms, report.error_residual_rms);
    fprintf(fid, 'Top absolute error pixels at best loop:\n');
    for k = 1:height(report.top_abs_error_pixels)
        fprintf(fid, '  row %2d col %2d: amp error %+9.6f, intensity %.6f\n', ...
            report.top_abs_error_pixels.row(k), ...
            report.top_abs_error_pixels.col(k), ...
            report.top_abs_error_pixels.amp_error(k), ...
            report.top_abs_error_pixels.intensity(k));
    end
    clear cleanup;

    makeAnalysisFigure(run_dir, diag_data, report, loops, best_raw, best_err, row_mean, col_mean, has_point_sampling, point_best_err);
    fprintf('Analysis report saved: %s\n', out_txt);
    fprintf('Analysis figure saved: %s\n', fullfile(run_dir, 'analysis_figure.png'));
end

function makeAnalysisFigure(run_dir, diag_data, report, loops, best_raw, best_err, row_mean, col_mean, has_point_sampling, point_best_err)
    fig = figure('Name', 'Hardware Diagnostics Analysis', ...
        'Color', 'w', ...
        'Visible', 'off', ...
        'Position', [80, 80, 1600, 1000]);
    tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(loops, diag_data.err_rms, 'o-', 'LineWidth', 1.4);
    hold on;
    if has_point_sampling && isfield(diag_data, 'point_err_rms')
        plot(loops, diag_data.point_err_rms, 's--', 'LineWidth', 1.1);
        legend({'area average', 'point sample'}, 'Location', 'northeast');
    end
    grid on;
    xlabel('Loop');
    ylabel('RMS');
    title(sprintf('Area RMS drop %.1f%%, best %02d', report.rms_drop_pct, report.best_loop));

    nexttile;
    yyaxis left;
    plot(loops, diag_data.sample_mean, 'o-', 'LineWidth', 1.2);
    ylabel('Mean intensity');
    yyaxis right;
    roi_std = squeeze(std(reshape(diag_data.sampled_roi_smooth_stack, [], numel(loops)), 0, 1));
    plot(loops, roi_std, 's-', 'LineWidth', 1.2);
    ylabel('Intensity std');
    grid on;
    xlabel('Loop');
    title('ROI intensity mean/std');

    nexttile;
    imagesc(best_raw);
    axis image;
    colorbar;
    title(sprintf('Best sampled ROI intensity, CV %.3f', report.best_intensity_cv));

    nexttile;
    if has_point_sampling
        imagesc(point_best_err);
        plot_title = sprintf('Point amp error, RMS %.3f', report.point_best_error_rms);
    else
        imagesc(best_err);
        plot_title = sprintf('Best amp error, p95 %.3f', report.best_error_abs_p95);
    end
    axis image;
    colorbar;
    title(plot_title);

    nexttile;
    plot(row_mean, 1:numel(row_mean), 'o-', 'LineWidth', 1.2);
    set(gca, 'YDir', 'reverse');
    grid on;
    xlabel('Mean amp error');
    ylabel('ROI row');
    title('Row-wise mean error');

    nexttile;
    plot(1:numel(col_mean), col_mean, 'o-', 'LineWidth', 1.2);
    grid on;
    xlabel('ROI col');
    ylabel('Mean amp error');
    title('Column-wise mean error');

    exportgraphics(fig, fullfile(run_dir, 'analysis_figure.png'), 'Resolution', 150);
    close(fig);
end

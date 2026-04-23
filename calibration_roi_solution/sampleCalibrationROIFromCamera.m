function [cam_roi, info] = sampleCalibrationROIFromCamera(img_raw, calib_file, zero_center_xy, row_start, row_end, col_start, col_end, algo_grid_size, area_sample_count)
% sampleCalibrationROIFromCamera samples the algorithm ROI directly from camera pixels.
%
% The algorithm ROI pixel centers are mapped through calib_data.affine_matrix
% into full-frame camera coordinates, then sampled with interp2. When
% area_sample_count > 1, each algorithm pixel is represented by the average
% of an area_sample_count-by-area_sample_count subgrid across that pixel cell.

    if nargin < 8 || isempty(algo_grid_size)
        algo_grid_size = 1080;
    end
    if nargin < 9 || isempty(area_sample_count)
        area_sample_count = 1;
    end
    area_sample_count = max(1, round(area_sample_count));

    if ~ismatrix(img_raw)
        error('img_raw must be a 2-D grayscale image.');
    end

    S = load(calib_file);
    if ~isfield(S, 'calib_data') || ~isfield(S.calib_data, 'affine_matrix')
        error('Calibration file must contain calib_data.affine_matrix.');
    end

    tform = makeAffineTransform(S.calib_data.affine_matrix);

    img_raw = double(img_raw);
    [cam_H, cam_W] = size(img_raw);

    [algo_col_grid, algo_row_grid] = meshgrid(col_start:col_end, row_start:row_end);
    algo_center = (algo_grid_size + 1) / 2;
    algo_x = algo_col_grid - algo_center;
    algo_y = algo_row_grid - algo_center;

    [center_dx, center_dy] = transformPointsForward(tform, algo_x, algo_y);
    center_sample_x = zero_center_xy(1) + center_dx;
    center_sample_y = zero_center_xy(2) + center_dy;

    if area_sample_count == 1
        sample_x = center_sample_x;
        sample_y = center_sample_y;
        in_bounds = sample_x >= 1 & sample_x <= cam_W & ...
                    sample_y >= 1 & sample_y <= cam_H;

        cam_roi = interp2(img_raw, sample_x, sample_y, 'linear', NaN);
        valid_samples = isfinite(cam_roi);
        cam_roi(~valid_samples) = 0;
        sample_x_range = [min(sample_x(:)), max(sample_x(:))];
        sample_y_range = [min(sample_y(:)), max(sample_y(:))];
        out_of_bounds_ratio = 1 - nnz(in_bounds) / numel(in_bounds);
        valid_pixel_ratio = nnz(valid_samples) / numel(valid_samples);
    else
        offsets = ((1:area_sample_count) - 0.5) / area_sample_count - 0.5;
        cam_sum = zeros(size(algo_x));
        valid_count = zeros(size(algo_x));
        in_bounds_count = zeros(size(algo_x));
        sample_x_min = inf;
        sample_x_max = -inf;
        sample_y_min = inf;
        sample_y_max = -inf;

        for oy = offsets
            for ox = offsets
                [cam_dx, cam_dy] = transformPointsForward(tform, algo_x + ox, algo_y + oy);
                sample_x = zero_center_xy(1) + cam_dx;
                sample_y = zero_center_xy(2) + cam_dy;
                in_bounds = sample_x >= 1 & sample_x <= cam_W & ...
                            sample_y >= 1 & sample_y <= cam_H;
                sampled = interp2(img_raw, sample_x, sample_y, 'linear', NaN);
                valid_samples = isfinite(sampled);
                sampled(~valid_samples) = 0;

                cam_sum = cam_sum + sampled;
                valid_count = valid_count + valid_samples;
                in_bounds_count = in_bounds_count + in_bounds;
                sample_x_min = min(sample_x_min, min(sample_x(:)));
                sample_x_max = max(sample_x_max, max(sample_x(:)));
                sample_y_min = min(sample_y_min, min(sample_y(:)));
                sample_y_max = max(sample_y_max, max(sample_y(:)));
            end
        end

        total_subsamples = area_sample_count^2;
        cam_roi = cam_sum ./ max(valid_count, 1);
        cam_roi(valid_count == 0) = 0;
        in_bounds = in_bounds_count == total_subsamples;
        valid_samples = valid_count == total_subsamples;
        sample_x_range = [sample_x_min, sample_x_max];
        sample_y_range = [sample_y_min, sample_y_max];
        out_of_bounds_ratio = 1 - sum(in_bounds_count(:)) / (numel(in_bounds_count) * total_subsamples);
        valid_pixel_ratio = sum(valid_count(:)) / (numel(valid_count) * total_subsamples);
        sample_x = center_sample_x;
        sample_y = center_sample_y;
    end

    info.out_of_bounds_ratio = out_of_bounds_ratio;
    info.valid_pixel_ratio = valid_pixel_ratio;
    info.sample_x_range = sample_x_range;
    info.sample_y_range = sample_y_range;
    info.sample_x = sample_x;
    info.sample_y = sample_y;
    info.in_bounds = in_bounds;
    info.all_subsamples_valid = valid_samples;
    info.area_sample_count = area_sample_count;
    info.area_subsample_count = area_sample_count^2;
    info.sampling_mode = 'point';
    if area_sample_count > 1
        info.sampling_mode = 'area_average';
    end
    info.algo_roi_size_hw = size(cam_roi);
    info.camera_size_hw = [cam_H, cam_W];
    info.mean_intensity = mean(cam_roi(:));
    info.total_intensity = sum(cam_roi(:));
end

function tform = makeAffineTransform(T)
    if isa(T, 'affine2d')
        tform = T;
        return;
    end

    if ~isequal(size(T), [3, 3])
        error('calib_data.affine_matrix must be an affine2d object or a 3x3 affine matrix.');
    end

    tform = affine2d(T);
end

function [cam_roi, info] = sampleCalibrationROIFromCamera(img_raw, calib_file, zero_center_xy, row_start, row_end, col_start, col_end, algo_grid_size)
% sampleCalibrationROIFromCamera samples the algorithm ROI directly from camera pixels.
%
% The algorithm ROI pixel centers are mapped through calib_data.affine_matrix
% into full-frame camera coordinates, then sampled with interp2. This keeps
% each feedback-error pixel aligned with the corresponding algorithm ROI pixel.

    if nargin < 8 || isempty(algo_grid_size)
        algo_grid_size = 1080;
    end

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

    [cam_dx, cam_dy] = transformPointsForward(tform, algo_x, algo_y);
    sample_x = zero_center_xy(1) + cam_dx;
    sample_y = zero_center_xy(2) + cam_dy;

    in_bounds = sample_x >= 1 & sample_x <= cam_W & ...
                sample_y >= 1 & sample_y <= cam_H;

    cam_roi = interp2(img_raw, sample_x, sample_y, 'linear', NaN);
    valid_samples = isfinite(cam_roi);
    cam_roi(~valid_samples) = 0;

    info.out_of_bounds_ratio = 1 - nnz(in_bounds) / numel(in_bounds);
    info.valid_pixel_ratio = nnz(valid_samples) / numel(valid_samples);
    info.sample_x_range = [min(sample_x(:)), max(sample_x(:))];
    info.sample_y_range = [min(sample_y(:)), max(sample_y(:))];
    info.sample_x = sample_x;
    info.sample_y = sample_y;
    info.in_bounds = in_bounds;
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

# MATLAB MCP Agent Log

- Task: Build a calibration-affine ROI version of the camera-in-the-loop flat-top beam script.
- Correction: The physical formula version is kept for comparison; this version uses the current `Calibration_Mapping_Data.mat` as the ROI authority.
- Required workflow: run `Pixel_Matching_Probe.m` before each experiment, then run `calibration_roi_solution/cam_in_loop.m`.
- Quality gates: max reprojection error <= 3 px, scale error <= 5%, rotation abs <= 5 deg.
- Non-hardware validation: run `simulate_calibration_flat_top_roi.m`.

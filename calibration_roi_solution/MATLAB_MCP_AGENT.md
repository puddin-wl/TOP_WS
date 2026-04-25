# MATLAB MCP Agent Log

- Task: Build a calibration-affine ROI version of the camera-in-the-loop flat-top beam script.
- Correction: The physical formula version is kept for comparison; this version uses the current `Calibration_Mapping_Data.mat` as the ROI authority.
- Required workflow: run `Pixel_Matching_Probe.m` before each experiment, then run `calibration_roi_solution/cam_in_loop.m`.
- Quality gates: max reprojection error <= 3 px, scale error <= 5%, rotation abs <= 5 deg.
- Non-hardware validation: run `simulate_calibration_flat_top_roi.m`.
- 2026-04-23 update: closed-loop error now uses calibration-affine reverse sampling via
  `sampleCalibrationROIFromCamera.m`, so the measured error matrix is pixel-aligned with the
  algorithm ROI instead of being produced by camera-box crop/rotate/resize.
- 2026-04-23 update: feedback sampling now supports area averaging. Hardware feedback uses
  a 7x7 sub-sample average per algorithm ROI pixel to better optimize camera-area
  uniformity across the physical ROI, while point sampling is saved as a diagnostic contrast.
- Working rule: after each code or diagnostic modification is completed and checked,
  choose an appropriate commit message, commit the change, and push the current branch.

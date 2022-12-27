# Accelerometer_Calibration_by_Kalman_Filtering

This task is about calibrating an accelerometer by designing a Kalman filter.

In this task, an object is assumed to move in one dimensional direction and the position and velocity of the object are estimated from model dynamics and biased measurements.

The main.m file is where the main program executed. In this program, one estimation source is provided and true model, IMU model and Kalman Filter parts are all constructed in this program to show the correctness of Kalman filtering. The errors of estimation of each state are graphed in this file.

In Orthogonality_Check.m and Orthogonalities_More.m files, you can check the correctness of Kalman filter implementation by examining the orthogonal properties of Kalman filter. 

In Sequtial_and_batch_update_KF.m file, it's demonstrated that sequential update and batch update are the same for Kalman filter when there are multiple measurement sources.

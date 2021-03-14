# systole_detection

T-wave end detection using Trapez area algorithm to detect the systolic window of the cardiac cycle (VÃ¡zquez-Seisdedos et al., 2011)

To segment the cardiac cycle into systole and diastole, one can compute the trial-specific phases based on cardio-mechanical events related to the ECG trace.

The ventricular systolic phase (further referred as "systole") is defined as a time between R peak of the QRS complex and the t-wave end, while diastole as the remaining part of the RR interval.

The trapez area algorithm is applied to encode the t-wave end in each trial. First, the t-peak is located as a local maximum within the physiologically plausible interval after the R peak containing the t-wave.

Subsequently, the algorithm computes a series of trapezes along the descending part of the t-wave signal, defining the point at which the trapeziumÂ's area gets maximal as the t-wave end.

Note: To use the code, change main_path variable to the path of systole_detection folder! 

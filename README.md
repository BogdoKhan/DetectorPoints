# DetectorPoints
A tool for the calculation of track registration point in multilayer (4 layer x 2 detectors) straw detector systems.
Input data: straw detector tube radius and radii of isochrones on registered detectors in form of vector 
(for example: {4.2, 0.0, 2.4, 0.0, 3.2, 0.0, 1.6, 0.0}),
where 0.0 or nan correspond to the detectors without signal.
Provides output in form of list of 3D registration points in the system of detectors.
Requires GSL.
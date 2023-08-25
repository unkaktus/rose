import numpy as np

# Create the transfer function as a series of triangles
def points_for_peaks(peaks, width, opacity):
    points = np.array([])
    for peak in peaks:
        new_points = [
            peak-width, 0.0, 0.5, 0.0,
            peak,  opacity, 0.5, 0.0,
            peak+width, 0.0, 0.5, 0.0,
        ]
        points = np.append(points, new_points)
    print(f'{opacity}')
    return points

peaks = np.arange(40, 111, 10)
print(peaks)

width = 0.1
opacity = 0.4
tf = GetOpacityTransferFunction('eSPL in dB')
tf.Points = points_for_peaks(peaks, width, opacity)
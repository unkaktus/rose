import numpy as np

def points_for_peaks(peaks, lwidth, opacity):
    points = np.array([])
    for peak in peaks:
        new_points = [
            peak-lwidth, 0.0, 0.5, 0.0,
            peak,  opacity, 0.5, 0.0,
            peak, 0.0, 0.5, 0.0,
        ]
        points = np.append(points, new_points)
    print(f'{opacity}')
    return points

peaks = np.arange(40, 111, 10)
print(peaks)
lwidth = 0.3
opacity = 0.7
tf = GetOpacityTransferFunction('eSPL in dB')
tf.Points = points_for_peaks(peaks, lwidth, opacity)
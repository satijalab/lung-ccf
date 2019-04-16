import os
import numpy as np
import pandas as pd
import re


class Skeleton:
    """
        For reading and analyzing skeleton analyses obtained
        with "Analyze Skeleton" plugin from Fiji.
    """
    def __init__(self, csv_file=None):
        assert os.path.isfile(csv_file)
        self.data_frame = pd.read_csv(csv_file)
        self.pre_process_skeleton()
        self.n_branches = len(self.data_frame)
        self.points = None
        self.set_skeleton_points()
        self.histogram = None
        self.set_histogram()

    def pre_process_skeleton(self):
        df = self.data_frame
        mask = df['Euclidean distance'] > 0
        df = df[mask]
        mask = df['running average length'] > 0
        df = df[mask]
        mask = df['average intensity'] == 255
        df = df[mask]
        mask = df['average intensity (inner 3rd)'] == 255
        df = df[mask]
        self.data_frame = df

    def set_skeleton_points(self):
        pts_1 = self.data_frame[['V1 x', 'V1 y', 'V1 z']].values
        pts_2 = self.data_frame[['V2 x', 'V2 y', 'V2 z']].values
        pts = np.vstack([pts_1, pts_2])
        self.points = pts

    def set_histogram(self):
        histogram = {}
        for point in self.points:
            pt_str = repr(point)
            if pt_str in histogram.keys():
                histogram[pt_str] += 1
            else:
                histogram[pt_str] = 1
        self.histogram = histogram

    def _generate_point_w_connectivity(self, connectivity):
        for pt_str, conn in self.histogram.items():
            if conn == connectivity:
                yield pt_str

    @staticmethod
    def _get_coordinates(pt_repr):
        coordinates_str = re.findall('\d+', pt_repr)
        coordinates = [float(coordinate_str) for coordinate_str in coordinates_str]
        return np.array(coordinates).astype('float')

    def get_points_w_connectivity(self, connectivity=1):
        list_data = list(self._generate_point_w_connectivity(connectivity))
        coordinates_list = [self._get_coordinates(pt) for pt in list_data]
        xs, ys, zs = np.stack(coordinates_list).T
        return xs, ys, zs



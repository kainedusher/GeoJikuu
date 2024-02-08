
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 4 16:11:45 2024

Title: hot_spot_analysis.py
Last Updated: GeoJikuu v0.27.53

Description:
This module contains classes for calculating autocorrelation metrics. 

    
Please refer to the official documentation for more information.

Author: Kaine Usher (kaine.usher1@gmail.com)
License: Apache 2.0 (see LICENSE file for details)

Copyright (c) 2024, Kaine Usher.
"""

import numpy as np
import pandas as pd
import math

# This version only uses binary weights
class GlobalMoranI:
    
    def __init__(self, data, coordinate_label):
        self.__data = data.to_dict()
        self.__coordinate_label = coordinate_label
        self.__results = {}
        
    def run(self, input_field, critical_distance=0, alpha=0.05, verbose=True):
        
        points = self.__data[self.__coordinate_label]
        
        # Calculate the total number of features
        n = len(points)
        
        # Compute the spatial weights based on critical distance
        n_panel = self.__find_neighbours(points, critical_distance, input_field)
        
        # Calculate the aggregate of all the spatial weights
        spatial_weights_agg = np.sum(n_panel["neighbours"])
        fraction_one = n / spatial_weights_agg
        
        # Multiply the columns element-wise and then sum the result
        df = pd.DataFrame(n_panel)
        fraction_two_numerator = (df['value_i_dev'] * df['value_j_dev'] * df['neighbours']).sum()
        
        # Find the sum of squared i value deviations
        unique_df = df.drop_duplicates(subset='i')
        fraction_two_denominator = ((unique_df['value_i_dev']) ** 2).sum()
        
        moran_i = fraction_one * fraction_two_numerator / fraction_two_denominator
        z_score = self.__compute_z(moran_i, n, n_panel, spatial_weights_agg)
        p_value = self.__p_value(z_score)
        
        return {"I": moran_i, "Z-SCORE": z_score, "P-VALUE": p_value}
        
        
    def __find_neighbours(self, points, critical_distance, input_field):
        
        point_values = self.__data[input_field]
        point_value_mean = np.mean(list(point_values.values()))
        
        n_panel = {"i": [], "j": [], "value_i_dev": [], "value_j_dev": [], "neighbours": []}
        
        for point_key_i in points:
            point_i = points[point_key_i]
            for point_key_j in points:
                point_j = points[point_key_j]
                
                n_panel["i"].append(point_i)
                n_panel["value_i_dev"].append(point_values[point_key_i] - point_value_mean)
                n_panel["j"].append(point_j)
                n_panel["value_j_dev"].append(point_values[point_key_j] - point_value_mean)

                if point_key_i == point_key_j:
                    n_panel["neighbours"].append(0)
                    continue
                
                distance_from_target = self.__euclidean_distance(point_i, point_j)
            
                if distance_from_target <= critical_distance:
                    n_panel["neighbours"].append(1)
                else:
                    n_panel["neighbours"].append(0)
        
        return n_panel
        
    def __compute_z(self, moran_i, n, n_panel, spatial_weights_agg):
        
        e_i = -1 / (n - 1)
        e_i_sqr = self.__compute_e_i_squared(n, n_panel, spatial_weights_agg)
        
        v_i = e_i_sqr - e_i**2
        
        return (moran_i - e_i) / (v_i)**0.5
    
    def __compute_e_i_squared(self, n, n_panel, s_0):
        
        df = pd.DataFrame(n_panel)
        
        # Compute s_1
        weight_sums = ((2* df['neighbours']) ** 2).sum() # This only works for binary 'symmetric' weights
        s_1 = 0.5 * weight_sums
        
        # Compute s_2
        matrix = self.__compute_weights_matrix(n_panel)
        s_2 = self.__compute_s2(np.array(matrix))
        
        # Compute D
        D = self.__compute_d(n, n_panel)
        
        # Compute A, B, C
        A = n * ((n**2 - 3*n + 3) * s_1 - n * s_2 + 3*s_0**2)
        B = D * ((n**2 - n) * s_1 - 2*n*s_2 + 6*s_0**2)
        C = (n-1) * (n-2) * (n-3) * s_0**2
        
        return (A - B) / C
    
    def __compute_weights_matrix(self, data):
                
        # Flatten the list of tuples and create a set of unique coordinates
        unique_coords = set(data['i']) | set(data['j'])  # Union of unique 'i' and 'j' coordinates
        
        # Create a mapping from coordinate to unique index
        coord_to_index = {coord: idx for idx, coord in enumerate(unique_coords)}
        
        # Determine the size of the adjacency matrix
        size = len(unique_coords)
        
        # Initialize the adjacency matrix with 0s (no edge)
        adj_matrix = [[0] * size for _ in range(size)]
        
        # Populate the adjacency matrix
        for i, j, adj in zip(data['i'], data['j'], data['neighbours']):
            # Convert coordinates to indexes and update the matrix based on adjacency
            i_index, j_index = coord_to_index[i], coord_to_index[j]
            adj_matrix[i_index][j_index] = 1 if adj else 0
            adj_matrix[j_index][i_index] = 1 if adj else 0  
            
        return adj_matrix
        
    def __compute_s2(self, W):
        n = W.shape[0]  # Assuming W is a square matrix
        total_sum = 0
        
        for i in range(n):
            inner_sum = 0
            for j in range(n):
                inner_sum += W[i, j] + W[j, i]
            total_sum += inner_sum ** 2
        
        return total_sum
    
    def __compute_d(self, n, n_panel):
        
        df = pd.DataFrame(n_panel)
        numerator = n * df['value_i_dev'].pow(4).sum()
        denominator = df['value_i_dev'].pow(2).sum()**2
        
        return numerator / denominator
        
    def __p_value(self, z_score):
        
        # Kudos to Sergei Winitzki
        # https://www.academia.edu/9730974/A_handy_approximation_for_the_error_function_and_its_inverse
        
        upper_bound = z_score * 10 / 2**0.5
        lower_bound = z_score / 2**0.5
        
        a = 8/(3*math.pi) * ((math.pi-3)/(4-math.pi))
        
        erf_upper = ((1 - math.exp(-upper_bound**2 * (4/math.pi+a*upper_bound**2) / (1+a*upper_bound**2)))**0.5)/2
        erf_lower = ((1 - math.exp(-lower_bound**2 * (4/math.pi+a*lower_bound**2) / (1+a*lower_bound**2)))**0.5)/2
        
        return 2 * (erf_upper - erf_lower)      
        
    
    def __euclidean_distance(self, x, y):

        if type(x) == str:
            x_string = x.strip('(').strip(')').split(", ")
            x = tuple([float(i) for i in x_string])

        if type(y) == str:
            y_string = y.strip('(').strip(')').split(", ")
            y = tuple([float(i) for i in y_string])
        
        euclid_distance = 0
    
        for i in range(0, len(x)):
            euclid_distance += (float(x[i]) - float(y[i]))**2
    
        return euclid_distance**0.5
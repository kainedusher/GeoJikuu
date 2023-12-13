# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:09:18 2023

Title: hot_spot_analysis.py
Last Updated: GeoJikuu v0.25.45

Description:
This module contains classes for performing hot spot analysis. 

    
Please refer to the official documentation for more information.

Author: Kaine Usher (kaine.usher1@gmail.com)
License: Apache 2.0 (see LICENSE file for details)

Copyright (c) 2023, Kaine Usher.
"""

import pandas as pd
import numpy as np
import math

class GiStarHotSpotAnalysis:
    
    def __init__(self, data, coordinate_label):
        self.__data = data.to_dict()
        self.__coordinate_label = coordinate_label
        self.__results = {}
        
    def run(self, input_field, critical_distance=0, alpha=0.05, verbose=True):
        
        results = {
            self.__coordinate_label: [],
            "neighbours": [],
            "z-score": [],
            "p-value": [],
            "significant": [],
            "type": []
            }
        
        points = self.__data[self.__coordinate_label]
        
        for target_key, target_coords in points.items():
            
            j_set = self.__find_neighbours(target_key, target_coords, points, critical_distance, input_field)
            z_score = self.__getis_ord_gi_star(j_set)
            p_value = self.__p_value(z_score)
            
            results[self.__coordinate_label].append(target_coords)
            results["neighbours"].append(j_set["neighbour"].count(True))
            results["z-score"].append(z_score)
            results["p-value"].append(p_value)
            
            if p_value*100 < alpha*100:
                results['significant'].append("TRUE")
            else:
                results['significant'].append("FALSE")
                
            if z_score >= 0:
                results['type'].append("HOT SPOT")
            else:
                results['type'].append("COLD SPOT")
                
        results = pd.DataFrame.from_dict(results)
        
        if verbose:
            significant_features = len(results[results['significant'] == "TRUE"])
            significant_hot = len(results[(results['significant'] == "TRUE") & (results['type'] == "HOT SPOT")])
            significant_cold = len(results[(results['significant'] == "TRUE") & (results['type'] == "COLD SPOT")])
            total_features = len(results)
            other_features = total_features - significant_features
            
            print("Getis-Ord Gi* Hot Spot Analysis Summary")
            print("---------------------------------------")
            print("Statistically Significant Features: " + str(significant_features))
            print("    Statistically Significant Hot Spots: " + str(significant_hot))
            print("    Statistically Significant Cold Spots: " + str(significant_cold))
            print("Non-Statistically Significant Features: " + str(other_features))
            print("Total Features: " + str(total_features))
                  
            print("")
            print("Null Hypothesis (H\N{SUBSCRIPT ZERO}): The observed pattern of the variable '" + str(input_field) + "' in feature \N{Double-Struck Italic Small I} is the result of spatial randomness alone.")
            print("Alpha Level (\N{GREEK SMALL LETTER ALPHA}): " + str(alpha))
            print("Critical Distance: " + str(critical_distance))
            print("Spatial Relationship Function: Inverse Distance")
            print("")
            
            if significant_features > 0:
                
                limit = 100 # This limits the number of results to be outputted in the sig results set.
                           # The number is arbitrary and might be changed in the future.
                sig_df = results[results['significant'] == "TRUE"]
                
                # Use slicing to limit the number of items
                limited_indices = [str(index) for index in sig_df.index][:limit]
                
                # Check if the total number of items exceeds the limit
                if len(sig_df) > limit:
                    limited_indices.append("...")
                
                # Join the items into a string
                sig_feature_labels_string = ', '.join(limited_indices)
                
                print("Verdict: Sufficient evidence to reject H\N{SUBSCRIPT ZERO} when \N{GREEK SMALL LETTER ALPHA} = " + str(alpha) + " for features \N{Double-Struck Italic Small I} \N{Element Of} {" + sig_feature_labels_string + "}")
            else:
                print("Verdict: Insufficient evidence to reject H\N{SUBSCRIPT ZERO} when \N{GREEK SMALL LETTER ALPHA} = " + str(alpha) + " for any of the analysed features.")
            
            
        return results
                                     
    
    def __getis_ord_gi_star(self, j_set):
        
        n = len(j_set["coords"])
        
        # Numerator
        sum_of_weights_times_values = 0
        for i in range(0, n):
            sum_of_weights_times_values += j_set["inverse_distance"][i] *  j_set["neighbour"][i] * j_set["value"][i]
        
        x_bar = sum(j_set["value"]) / n
        
        sum_of_weights = 0
        for i in range(0, n):
            sum_of_weights += j_set["inverse_distance"][i] *  j_set["neighbour"][i]
        
        numerator = sum_of_weights_times_values - x_bar * sum_of_weights
        
        # Denominator
        sum_of_squared_values = 0
        for i in range(0, n):
            sum_of_squared_values += j_set["value"][i]**2
            
        s = ((sum_of_squared_values / n) - (x_bar)**2)**0.5
        
        sum_of_squared_weights = 0
        for i in range(0, n):
            sum_of_squared_weights += (j_set["inverse_distance"][i] *  j_set["neighbour"][i])**2
        
        denominator = s * ((n * sum_of_squared_weights - (sum_of_weights)**2)/(n - 1))**0.5
        
        gets_ord_gi_star = numerator / denominator
        
        return gets_ord_gi_star
                
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
        
    
    def __find_neighbours(self, target_key, target_coords, points, critical_distance, input_field):
        
        point_values = self.__data[input_field]
        j_set = {"coords": [], "distance": [], "inverse_distance": [], "value": [], "neighbour": []}
        
        index = 0
        for point_key, point_coords in points.items():
            
            distance_from_target = self.__euclidean_distance(target_coords, point_coords)
            j_set["coords"].append(point_coords)
            j_set["distance"].append(distance_from_target)
            
            if distance_from_target == 0:
                j_set["inverse_distance"].append(1)
            else:
                j_set["inverse_distance"].append(1 / distance_from_target)
                
            j_set["value"].append(point_values[index])
            
            if distance_from_target <= critical_distance:
                j_set["neighbour"].append(True)
            else:
                j_set["neighbour"].append(False)
                
            index += 1
        
        
        return j_set
                
    
class STGiStarHotSpotAnalysis:
    
    def __init__(self, data, coordinate_label, time_label=None):
        self.__data = data.to_dict()
        
        if time_label is not None:
            self.__coordinate_label = self.__combine_labels(coordinate_label, time_label)
        else:
            self.__coordinate_label = coordinate_label
        self.__results = {}
        
    def run(self, input_field, critical_distance=0, critical_time=None, alpha=0.05, verbose=True):
        
        results = {
            self.__coordinate_label: [],
            "neighbours": [],
            "z-score": [],
            "p-value": [],
            "significant": [],
            "type": []
            }
        
        points = self.__data[self.__coordinate_label]
        
        for target_key, target_coords in points.items():
            
            j_set = self.__find_neighbours(target_key, target_coords, points, critical_distance, critical_time, 
                                           input_field)
            z_score = self.__getis_ord_gi_star(j_set)
            p_value = self.__p_value(z_score)
            
            results[self.__coordinate_label].append(target_coords)
            results["neighbours"].append(j_set["neighbour"].count(True))
            results["z-score"].append(z_score)
            results["p-value"].append(p_value)
            
            if p_value*100 < alpha*100:
                results['significant'].append("TRUE")
            else:
                results['significant'].append("FALSE")
                
            if z_score >= 0:
                results['type'].append("HOT SPOT")
            else:
                results['type'].append("COLD SPOT")
                
        results = pd.DataFrame.from_dict(results)
            
        if verbose:
            significant_features = len(results[results['significant'] == "TRUE"])
            significant_hot = len(results[(results['significant'] == "TRUE") & (results['type'] == "HOT SPOT")])
            significant_cold = len(results[(results['significant'] == "TRUE") & (results['type'] == "COLD SPOT")])
            total_features = len(results)
            other_features = total_features - significant_features
            
            print("Getis-Ord Gi* Hot Spot Analysis Summary")
            print("---------------------------------------")
            print("Statistically Significant Features: " + str(significant_features))
            print("    Statistically Significant Hot Spots: " + str(significant_hot))
            print("    Statistically Significant Cold Spots: " + str(significant_cold))
            print("Non-Statistically Significant Features: " + str(other_features))
            print("Total Features: " + str(total_features))
                  
            print("")
            print("Null Hypothesis (H\N{SUBSCRIPT ZERO}): The observed pattern of the variable '" + str(input_field) + "' in feature \N{Double-Struck Italic Small I} is the result of spatiotemporal randomness alone.")
            print("Alpha Level (\N{GREEK SMALL LETTER ALPHA}): " + str(alpha))
            print("Critical Distance: " + str(critical_distance))
            print("Spatial Relationship Function: Inverse Spacetime Distance")
            print("")
            
            if significant_features > 0:
                
                limit = 100 # This limits the number of results to be outputted in the sig results set.
                           # The number is arbitrary and might be changed in the future.
                sig_df = results[results['significant'] == "TRUE"]
                
                # Use slicing to limit the number of items
                limited_indices = [str(index) for index in sig_df.index][:limit]
                
                # Check if the total number of items exceeds the limit
                if len(sig_df) > limit:
                    limited_indices.append("...")
                
                # Join the items into a string
                sig_feature_labels_string = ', '.join(limited_indices)
                
                print("Verdict: Sufficient evidence to reject H\N{SUBSCRIPT ZERO} when \N{GREEK SMALL LETTER ALPHA} = " + str(alpha) + " for features \N{Double-Struck Italic Small I} \N{Element Of} {" + sig_feature_labels_string + "}")
            else:
                print("Verdict: Insufficient evidence to reject H\N{SUBSCRIPT ZERO} when \N{GREEK SMALL LETTER ALPHA} = " + str(alpha) + " for any of the analysed features.")
            
            
        return results
    
    def __getis_ord_gi_star(self, j_set):
        
        n = len(j_set["coords"])
        
        # Numerator
        sum_of_weights_times_values = 0
        for i in range(0, n):
            sum_of_weights_times_values += j_set["inverse_spacetime_distance"][i] *  j_set["neighbour"][i] * j_set["value"][i]
        
        x_bar = sum(j_set["value"]) / n
        
        sum_of_weights = 0
        for i in range(0, n):
            sum_of_weights += j_set["inverse_spacetime_distance"][i] *  j_set["neighbour"][i]
        
        numerator = sum_of_weights_times_values - x_bar * sum_of_weights
        
        # Denominator
        sum_of_squared_values = 0
        for i in range(0, n):
            sum_of_squared_values += j_set["value"][i]**2
            
        s = ((sum_of_squared_values / n) - (x_bar)**2)**0.5
        
        sum_of_squared_weights = 0
        for i in range(0, n):
            sum_of_squared_weights += (j_set["inverse_spacetime_distance"][i] *  j_set["neighbour"][i])**2
        
        denominator = s * ((n * sum_of_squared_weights - (sum_of_weights)**2)/(n - 1))**0.5
        
        gets_ord_gi_star = numerator / denominator
        
        return gets_ord_gi_star
    
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
    
    def __find_neighbours(self, target_key, target_coords, points, critical_distance, critical_time, input_field):
        
        point_values = self.__data[input_field]
        j_set = {"coords": [], "spacetime_distance": [], "inverse_spacetime_distance": [], "value": [], "neighbour": []}
        
        index = 0
        for point_key, point_coords in points.items():
            
            spatial_distance_from_target = self.__euclidean_distance(target_coords[:-1], point_coords[:-1])
            temporal_distance_from_target = abs(target_coords[len(target_coords)-1] - point_coords[len(point_coords)-1])
            spacetime_distance_from_target = self.__euclidean_distance(target_coords, point_coords)
            
            j_set["coords"].append(point_coords)
            j_set["spacetime_distance"].append(spacetime_distance_from_target)
            
            if spacetime_distance_from_target == 0:
                j_set["inverse_spacetime_distance"].append(1)
            else:
                j_set["inverse_spacetime_distance"].append(1 / spacetime_distance_from_target)
                
            j_set["value"].append(point_values[index])
            
            if critical_time is not None:
                if spatial_distance_from_target <= critical_distance and temporal_distance_from_target <= critical_time:
                    j_set["neighbour"].append(True)
                else:
                    j_set["neighbour"].append(False)
            else:
                if spatial_distance_from_target <= critical_distance:
                    j_set["neighbour"].append(True)
                else:
                    j_set["neighbour"].append(False)
                
            index += 1
        
        
        return j_set
    
    def __combine_labels(self, coordinate_label, time_label):
        
        spatial_coordinates = self.__data[coordinate_label]
        temporal_coordinates = self.__data[time_label]
        spatial_temporal_coordinates = []
        
        for i in range(0, len(spatial_coordinates)):
            tuple_list = []
            for coord in spatial_coordinates[i]:
                tuple_list.append(coord)
            tuple_list.append(temporal_coordinates[i])
            spatial_temporal_coordinates.append(i)
            spatial_temporal_coordinates.append(tuple(tuple_list))

        res_dct = {spatial_temporal_coordinates[i]: spatial_temporal_coordinates[i + 1] for i in range(0, len(spatial_temporal_coordinates), 2)}
        self.__data["st_coordinates"] = res_dct
        return "st_coordinates"
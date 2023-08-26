# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:09:18 2023

Title: hot_spot_analysis.py
Last Updated: GeoJikuu v0.23.31

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
        
        
    def run(self, input_field, alpha=0.05, verbose=True):
        
        results = {
            self.__coordinate_label: [],
            "z-score": [],
            "p-value": [],
            "significant": [],
            "type": []
            }
        
        j_set = self.__data[input_field]
        points = self.__data[self.__coordinate_label]
        
        for key, value in points.items():
            
            z_score = self.__getis_ord_gi_star(value, j_set, points)
            p_value = self.__p_value(z_score)
            
            results[self.__coordinate_label].append(value)
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
            significant_clusters = len(results[results['significant'] == "TRUE"])
            significant_hot = len(results[(results['significant'] == "TRUE") & (results['type'] == "HOT SPOT")])
            significant_cold = len(results[(results['significant'] == "TRUE") & (results['type'] == "COLD SPOT")])
            total_clusters = len(results)
            other_clusters = total_clusters - significant_clusters
            
            print("Getis-Ord Gi* Hot Spot Analysis Summary")
            print("---------------------------------------")
            print("Statistically Significant Clusters: " + str(significant_clusters))
            print("    Statistically Significant Hot Spots: " + str(significant_hot))
            print("    Statistically Significant Cold Spots: " + str(significant_cold))
            print("Non-Statistically Significant Clusters: " + str(other_clusters))
            print("Total Clusters: " + str(total_clusters))
                  
            print("")
            print("Null Hypothesis (H\N{SUBSCRIPT ZERO}): The observed pattern of the variable '" + str(input_field) + "' in cluster \N{Double-Struck Italic Small I} is the result of spatial randomness alone.")
            print("Alpha Level (\N{GREEK SMALL LETTER ALPHA}): " + str(alpha))
            print("")
            
            if significant_clusters > 0:
                sig_df = results[results['significant'] == "TRUE"]
                sig_cluster_labels_string = ', '.join([str(index) for index in sig_df.index])
                print("Verdict: Sufficient evidence to reject H\N{SUBSCRIPT ZERO} when \N{GREEK SMALL LETTER ALPHA} = " + str(alpha) + " for clusters \N{Double-Struck Italic Small I} = {" + sig_cluster_labels_string + "}")
            else:
                print("Verdict: Insufficient evidence to reject H\N{SUBSCRIPT ZERO} when \N{GREEK SMALL LETTER ALPHA} = " + str(alpha) + " for any of the analysed clusters.")
            
            
        return results
                                      
    def __getis_ord_gi_star(self, i_point, j_set, points):
        
        x_bar = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            x_bar = x_bar + j_set[i]
       
        x_bar = x_bar / len(j_set)
        
        s = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            s = s + j_set[i]**2
        
        s = ((s / len(j_set)) - x_bar**2)**0.5
        
        gi_star_num_sum_one = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            gi_star_num_sum_one = gi_star_num_sum_one + (self.__euclidean_distance(i_point, points[i]) * j_set[i])
        
        gi_star_num_sum_two = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            gi_star_num_sum_two = gi_star_num_sum_two + (self.__euclidean_distance(i_point, points[i]))
            
        gi_star_num = gi_star_num_sum_one - x_bar * gi_star_num_sum_two
        
        gi_star_den_sum_one = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            gi_star_den_sum_one = gi_star_den_sum_one + (self.__euclidean_distance(i_point, points[i])**2)
            
        gi_star_den_sum_two = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            gi_star_den_sum_two = gi_star_den_sum_two + (self.__euclidean_distance(i_point, points[i]))
            
        gi_star_den = s * ((len(j_set) * gi_star_den_sum_one - (gi_star_den_sum_two)**2) /(len(j_set) - 1))**0.5
        
        gi_star = gi_star_num / gi_star_den
        
        return gi_star
        
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
        
class STGiStarHotSpotAnalysis:
    
    def __init__(self, data, coordinate_label):
        self.__data = data.to_dict()
        self.__coordinate_label = coordinate_label
        self.__results = {}
        
    def run(self, input_field, alpha=0.05, verbose=True):
        
        results = {
            self.__coordinate_label: [],
            "z-score": [],
            "p-value": [],
            "significant": [],
            "type": []
            }
        
        j_set = self.__data[input_field]
        points = self.__data[self.__coordinate_label]
        
        for key, value in points.items():
            
            z_score = self.__getis_ord_gi_star(value, j_set, points)
            p_value = self.__p_value(z_score)
            
            results[self.__coordinate_label].append(value)
            results["z-score"].append(z_score)
            results["p-value"].append(p_value)
            
            if p_value*100 < alpha*100:
                results["significant"].append("TRUE")
            else:
                results["significant"].append("FALSE")
                
            if z_score >= 0:
                results["type"].append("HOT SPOT")
            else:
                results["type"].append("COLD SPOT")
                
        
        results = pd.DataFrame.from_dict(results)
        
        if verbose:
            significant_clusters = len(results[results['significant'] == "TRUE"])
            significant_hot = len(results[(results['significant'] == "TRUE") & (results['type'] == "HOT SPOT")])
            significant_cold = len(results[(results['significant'] == "TRUE") & (results['type'] == "COLD SPOT")])
            total_clusters = len(results)
            other_clusters = total_clusters - significant_clusters
            
            print("Spacetime Getis-Ord Gi* Hot Spot Analysis Summary")
            print("-------------------------------------------------")
            print("Statistically Significant Clusters: " + str(significant_clusters))
            print("    Statistically Significant Hot Spots: " + str(significant_hot))
            print("    Statistically Significant Cold Spots: " + str(significant_cold))
            print("Non-Statistically Significant Clusters: " + str(other_clusters))
            print("Total Clusters: " + str(total_clusters))
                  
            print("")
            print("Null Hypothesis (H\N{SUBSCRIPT ZERO}): The observed pattern of the variable '" + str(input_field) + "' in cluster \N{Double-Struck Italic Small I} is the result of spatiotemporal randomness alone.")
            print("Alpha Level (\N{GREEK SMALL LETTER ALPHA}): " + str(alpha))
            print("")
            
            if significant_clusters > 0:
                sig_df = results[results['significant'] == "TRUE"]
                sig_cluster_labels_string = ', '.join([str(index) for index in sig_df.index])
                print("Verdict: Sufficient evidence to reject H\N{SUBSCRIPT ZERO} when \N{GREEK SMALL LETTER ALPHA} = " + str(alpha) + " for clusters \N{Double-Struck Italic Small I} = {" + sig_cluster_labels_string + "}")
            else:
                print("Verdict: Insufficient evidence to reject H\N{SUBSCRIPT ZERO} when \N{GREEK SMALL LETTER ALPHA} = " + str(alpha) + " for any of the analysed clusters.")
                 
        return results
    
    def __getis_ord_gi_star(self, i_point, j_set, points):
        
        x_bar = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            x_bar = x_bar + j_set[i]
       
        x_bar = x_bar / len(j_set)
        
        s = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            s = s + j_set[i]**2
        
        s = ((s / len(j_set)) - x_bar**2)**0.5
        
        gi_star_num_sum_one = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            gi_star_num_sum_one = gi_star_num_sum_one + (self.__euclidean_distance(i_point, points[i]) * j_set[i])
        
        gi_star_num_sum_two = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            gi_star_num_sum_two = gi_star_num_sum_two + (self.__euclidean_distance(i_point, points[i]))
            
        gi_star_num = gi_star_num_sum_one - x_bar * gi_star_num_sum_two
        
        gi_star_den_sum_one = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            gi_star_den_sum_one = gi_star_den_sum_one + (self.__euclidean_distance(i_point, points[i])**2)
            
        gi_star_den_sum_two = 0
        
        for i in range(0, len(j_set)):
            if points[i] == i_point:
                continue
            gi_star_den_sum_two = gi_star_den_sum_two + (self.__euclidean_distance(i_point, points[i]))
            
        gi_star_den = s * ((len(j_set) * gi_star_den_sum_one - (gi_star_den_sum_two)**2) /(len(j_set) - 1))**0.5
        
        gi_star = gi_star_num / gi_star_den
        
        return gi_star
    
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
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 15:09:18 2023

Title: hot_spot_analysis.py
Last Updated: GeoJikuu vx.x.x

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
        
        
    def run(self, input_field, alpha=0.05):
        
        results = {
            self.__coordinate_label: [],
            "Z-score": [],
            "p-value": [],
            "Significant": [],
            "Type": []
            }
        
        j_set = self.__data[input_field]
        points = self.__data[self.__coordinate_label]
        
        for key, value in points.items():
            
            z_score = self.__getis_ord_gi_star(value, j_set, points)
            p_value = self.__p_value(z_score)
            
            results[self.__coordinate_label].append(value)
            results["Z-score"].append(z_score)
            results["p-value"].append(p_value)
            
            if p_value*100 < alpha*100:
                results["Significant"].append("True")
            else:
                results["Significant"].append("False")
                
            if z_score >= 0:
                results["Type"].append("Hot Spot")
            else:
                results["Type"].append("Cold Spot")
                
        self.__results = results
                 
        return pd.DataFrame.from_dict(results)
                                      
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
        
    def run(self, input_field, alpha=0.05):
        
        results = {
            self.__coordinate_label: [],
            "Z-score": [],
            "p-value": [],
            "Significant": [],
            "Type": []
            }
        
        j_set = self.__data[input_field]
        points = self.__data[self.__coordinate_label]
        
        for key, value in points.items():
            
            z_score = self.__getis_ord_gi_star(value, j_set, points)
            p_value = self.__p_value(z_score)
            
            results[self.__coordinate_label].append(value)
            results["Z-score"].append(z_score)
            results["p-value"].append(p_value)
            
            if p_value*100 < alpha*100:
                results["Significant"].append("True")
            else:
                results["Significant"].append("False")
                
            if z_score >= 0:
                results["Type"].append("Hot Spot")
            else:
                results["Type"].append("Cold Spot")
                
        self.__results = results
                 
        return pd.DataFrame.from_dict(results)
    
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
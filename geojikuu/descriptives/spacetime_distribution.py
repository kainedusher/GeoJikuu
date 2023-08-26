# -*- coding: utf-8 -*-
"""
Created on Tue May 16 20:18:34 2023
Last Updated: GeoJikuu v0.23.31

Description:
This module contains classes for examining the shape, centre, and spread of 
spatiotemporal data.

All distance results are given in metres (or metres squared, in the case of variance).
Temporal output is given in days or in date format dd/mm/yyyy.
    
Please refer to the official documentation for more information.

Author: Kaine Usher (kaine.usher1@gmail.com)
License: Apache 2.0 (see LICENSE file for details)

Copyright (c) 2023, Kaine Usher.
"""

import math
import numpy as np
from datetime import *

class SpacetimePointDistribution:
    
    def __init__(self, lat_long_time_coordinates, date_format = "%d/%m/%Y"): 
        self.lat_long_time_coordinates = lat_long_time_coordinates
        self.date_format = date_format
        
    def count_points(self):
        return len(self.lat_long_time_coordinates)
    
    def mean_displacement(self, reference_point = ()):
        
        total_distances = 0
        total_time_differences = 0
        n = self.count_points()
        
        pair_count = 0
        if reference_point == (): 
            for i, point1 in enumerate(self.lat_long_time_coordinates):
                for j, point2 in enumerate(self.lat_long_time_coordinates):
                    if i != j and i < j:
                        spatial_point1 = (point1[0], point1[1])
                        spatial_point2 = (point2[0], point2[1])
                        total_distances += self.__haversine(spatial_point1, spatial_point2)
                        total_time_differences += self.__temporal_difference(point1[2], point2[2])
                        pair_count = pair_count + 1
                        
            mean_displacement_dict = {"MEAN SPATIAL DISPLACEMENT": (total_distances / pair_count), 
                                      "MEAN TEMPORAL DISPLACEMENT": (total_time_differences / pair_count)}
            return mean_displacement_dict
                        
        for lat_long_time_coordinate in self.lat_long_time_coordinates:
            spatial_point = (lat_long_time_coordinate[0], lat_long_time_coordinate[1])
            spatial_point_reference = (reference_point[0], reference_point[1])
            total_distances += self.__haversine(spatial_point_reference, spatial_point)
            
            total_time_differences += self.__temporal_difference(reference_point[2], lat_long_time_coordinate[2])
            mean_displacement_dict = {"MEAN SPATIAL DISPLACEMENT": (total_distances / n),
                                      "MEAN TEMPORAL DISPLACEMENT": (total_time_differences / n)}
            
        return mean_displacement_dict
            
    
    def displacement_std(self, reference_point = (), population=False):
        
        distances = []
        temporal_differences = []
        displacement_mean = 0
        temporal_difference_mean = 0
        n = 0
        
        if reference_point == (): 
            mean_displacement_dict = self.mean_displacement()
            displacement_mean = mean_displacement_dict["MEAN SPATIAL DISPLACEMENT"]
            temporal_difference_mean = mean_displacement_dict["MEAN TEMPORAL DISPLACEMENT"]
            
            for i, point1 in enumerate(self.lat_long_time_coordinates):
                for j, point2 in enumerate(self.lat_long_time_coordinates):
                    if i != j and i < j:
                        spatial_point1 = (point1[0], point1[1])
                        spatial_point2 = (point2[0], point2[1])
                        distances.append(self.__haversine(spatial_point1, spatial_point2)) 
                        temporal_differences.append(self.__temporal_difference(point1[2], point2[2]))
                        n = n + 1
        else:
            mean_displacement_dict = self.mean_displacement(reference_point)
            displacement_mean = mean_displacement_dict["MEAN SPATIAL DISPLACEMENT"]
            temporal_difference_mean = mean_displacement_dict["MEAN TEMPORAL DISPLACEMENT"]
            
            n = self.count_points()
            for lat_long_time_coordinate in self.lat_long_time_coordinates:
                distances.append(self.__haversine((reference_point[0], reference_point[1]), (lat_long_time_coordinate[0], lat_long_time_coordinate[1])))
                temporal_differences.append(self.__temporal_difference(reference_point[2], lat_long_time_coordinate[2]))
        
        residuals_sum_dist = 0
        residuals_sum_temp = 0
        for x in distances:
            residuals_sum_dist += (x - displacement_mean)**2
        
        for x in temporal_differences:
            residuals_sum_temp += (x - temporal_difference_mean)**2
        
        if population:
            dist_std = (residuals_sum_dist / n)**0.5
            temp_std = (residuals_sum_temp / n)**0.5
        else:
            dist_std = (residuals_sum_dist / (n-1))**0.5
            temp_std = (residuals_sum_temp / (n-1))**0.5
        
        displacement_std_dict = {"SPATIAL DISPLACEMENT STD": dist_std,
                                 "TEMPORAL DISPLACEMENT STD": temp_std}
        
        return displacement_std_dict
    
    
    def displacement_variance(self, reference_point = (), population=False):
        
        
        if population:
            displacement_std_dict = self.displacement_std(reference_point, population=True)
        else:
            displacement_std_dict = self.displacement_std(reference_point)
    
        displacement_variance_dict = {"SPATIAL DISPLACEMENT VARIANCE": displacement_std_dict["SPATIAL DISPLACEMENT STD"]**2,
                                      "TEMPORAL DISPLACEMENT VARIANCE": displacement_std_dict["TEMPORAL DISPLACEMENT STD"]**2}
        
        return displacement_variance_dict
    
    def displacement_quartiles(self, reference_point=()):
        
        distances = []
        temporal_differences = []
        
        if reference_point == (): 
            for i, point1 in enumerate(self.lat_long_time_coordinates):
                for j, point2 in enumerate(self.lat_long_time_coordinates):
                    if i != j and i < j:
                        spatial_point1 = (point1[0], point1[1])
                        spatial_point2 = (point2[0], point2[1])
                        distances.append(self.__haversine(spatial_point1, spatial_point2)) 
                        temporal_differences.append(self.__temporal_difference(point1[2], point2[2]))
        else:    
            for lat_long_time_coordinate in self.lat_long_time_coordinates:
                distances.append(self.__haversine((reference_point[0], reference_point[1]), (lat_long_time_coordinate[0], lat_long_time_coordinate[1])))
                temporal_differences.append(self.__temporal_difference(reference_point[2], lat_long_time_coordinate[2]))
                
        dis_minimum = np.min(distances)
        dis_q1 = np.quantile(distances, 0.25)
        dis_median = np.median(distances)
        dis_q3 = np.quantile(distances, 0.75)
        dis_maximum = np.max(distances)
        dis_iqr = dis_q3 - dis_q1
        dis_range = dis_maximum - dis_minimum
        
        temp_minimum = np.min(temporal_differences)
        temp_q1 = np.quantile(temporal_differences, 0.25)
        temp_median = np.median(temporal_differences)
        temp_q3 = np.quantile(temporal_differences, 0.75)
        temp_maximum = np.max(temporal_differences)
        temp_iqr = temp_q3 - temp_q1
        temp_range = temp_maximum - temp_minimum
        
        return {"MIN": {"DIS": dis_minimum, "TEMP": temp_minimum}, "Q1": {"DIS": dis_q1, "TEMP": temp_q1}, 
                "MEDIAN": {"DIS": dis_median, "TEMP": temp_median}, "Q3": {"DIS": dis_q3, "TEMP": temp_q3}, 
                "MAX": {"DIS": dis_maximum, "TEMP": temp_maximum}, "IQR": {"DIS": dis_iqr, "TEMP": temp_iqr}, 
                "RANGE": {"DIS": dis_range, "TEMP": temp_range}}
    

    def geo_temporal_midpoint(self):
        x_coords = []
        y_coords = []
        z_coords = []
        
        days = []
        
        for lat_long_time_coordinate in self.lat_long_time_coordinates:
            cartesian_coord = self.__latlong_to_cartesian((lat_long_time_coordinate[0], lat_long_time_coordinate[1]))
            x_coords.append(cartesian_coord[0])
            y_coords.append(cartesian_coord[1])
            z_coords.append(cartesian_coord[2])
            days.append(self.__date_to_days(lat_long_time_coordinate[2]))
            
        geo_mean = self.__cartesian_to_latlong((np.mean(x_coords), np.mean(y_coords), np.mean(z_coords)))
        temporal_mean = self.__days_to_date(np.mean(days))
        return ((geo_mean[0], geo_mean[1], temporal_mean))
        

    # Returns the difference in days
    def __temporal_difference(self, time1, time2):
        
        f_date = datetime.strptime(time2, self.date_format).date()
        l_date = datetime.strptime(time1, self.date_format).date()
        delta = l_date - f_date     
        
        return (abs(delta.days))
    
    def __date_to_days(self, date):
        
        return self.__temporal_difference(date, datetime.strptime("01/01/0001", "%d/%m/%Y").strftime(self.date_format))
    
    
    def __days_to_date(self, days):
        
        return (datetime.strptime("01/01/0001", "%d/%m/%Y").date() + timedelta(days=days)).strftime(self.date_format)
    
    # Input: An angle in degrees
    # Output: The input angle in radians
    def __degrees_to_rads(self, value):
        return (value * math.pi)/180
    
    # Input: An angle in radians
    # Output: The input angle in degrees
    def __rads_to_degrees(self, value):
        return (value * 180)/math.pi
    
    # Input: Two points in the form (Lat, Long)
    # Output: The displacement between those two points in m.
    def __haversine(self, p1, p2):
        
        p1_lat = self.__degrees_to_rads(p1[0])
        #p1_long = self.__degrees_to_rads(p1[1])
        p2_lat = self.__degrees_to_rads(p2[0])
        #p2_long = self.__degrees_to_rads(p2[1])
        p2_p1_lat_delta = self.__degrees_to_rads(p2[0] - p1[0])
        p2_p1_long_delta = self.__degrees_to_rads(p2[1] - p1[1])
        
        a = math.pow(math.sin((p2_p1_lat_delta)/2), 2) + math.cos(p1_lat) * math.cos(p2_lat) * math.pow(math.sin((p2_p1_long_delta)/2), 2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        return 6371 * c
    
    # Input: A point in the form (Lat, Long)
    # Output: The Cartesian equivalent to the input point; in the form (x, y, z)
    def __latlong_to_cartesian(self, lat_long_coordinate = ()):
        lat = self.__degrees_to_rads(lat_long_coordinate[0])
        long = self.__degrees_to_rads(lat_long_coordinate[1])
        
        x = math.cos(lat) * math.cos(long)
        y = math.cos(lat) * math.sin(long)
        z = math.sin(lat)
        
        return (x, y, z) 
    
    # Input: A Cartesian point in the form (x, y, z)
    # Output: The (Lat, Long) equivalent to the input point
    def __cartesian_to_latlong(self, cartesian_coordinate = ()):
        x = cartesian_coordinate[0]
        y = cartesian_coordinate[1]
        z = cartesian_coordinate[2]
        
        long = math.atan2(y, x)
        hyp = math.sqrt(x * x + y * y)
        lat = math.atan2(z, hyp)
        
        return (self.__rads_to_degrees(lat), self.__rads_to_degrees(long))
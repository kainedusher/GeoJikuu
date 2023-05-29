# -*- coding: utf-8 -*-
"""
Created on Tue May  9 14:11:40 2023
Last Updated: GeoJikuu v0.1.0

Description:
This module contains classes for examining the shape, centre, and spread of 
geospatial data.

All results are given in metres (or metres squared, in the case of variance). 
    
Please refer to the official documentation for more information.

Author: Kaine Usher (kaine.usher1@gmail.com)
License: Apache 2.0 (see LICENSE file for details)

Copyright (c) 2023, Kaine Usher.
"""
import math
import numpy as np

class PointDistribution:
    # Input: A list of (Lat, Long) coordinates
    # Output: N/A. This is the constructor
    def __init__(self, lat_long_coordinates): 
        self.lat_long_coordinates = lat_long_coordinates
        
    # Input: None
    # Output: The number of points in the dataset
    def count_points(self):
        return len(self.lat_long_coordinates)
    
    # Input: An optional reference point of the form (lat, long)
    # Output: The mean displacement from the reference point to all other points. If no
    # reference point is provided, then the mean displacment of each pair of points
    # is calculated. 
    def mean_displacement(self, reference_point = ()):
        
        total_distances = 0
        n = self.count_points()
        
        pair_count = 0
        if reference_point == (): 
            for i, point1 in enumerate(self.lat_long_coordinates):
                for j, point2 in enumerate(self.lat_long_coordinates):
                    if i != j and i < j:
                        total_distances += self.__haversine(point1, point2)
                        pair_count = pair_count + 1
            return total_distances / pair_count
                        
        for lat_long_coordinate in self.lat_long_coordinates:
            total_distances += self.__haversine(reference_point, lat_long_coordinate)
            
        return (total_distances / n)
    
    # Input: A point of the form (lat, long)
    # Output: The standard deviation of the displacements from the input point 
    #         to all other points. Provi
    def displacement_std(self, reference_point = (), population=False):
        
        distances = []
        mean = 0
        n = 0
        
        if reference_point == (): 
            mean = self.mean_displacement()
            for i, point1 in enumerate(self.lat_long_coordinates):
                for j, point2 in enumerate(self.lat_long_coordinates):
                    if i != j and i < j:
                        distances.append(self.__haversine(point1, point2)) 
                        n = n + 1
        else:    
            mean = self.mean_displacement(reference_point)
            n = self.count_points()
            for lat_long_coordinate in self.lat_long_coordinates:
                distances.append(self.__haversine(reference_point, lat_long_coordinate))
        
        residuals_sum = 0
        for x in distances:
            residuals_sum = residuals_sum + (x - mean)**2
        
        if population:
            std = (residuals_sum / n)**0.5
        else:
            std = (residuals_sum / (n-1))**0.5
        
        return std
    
    def displacement_variance(self, reference_point = (), population=False):
        
        if population:
            return self.displacement_std(reference_point, population=True) ** 2
        
        return self.displacement_std(reference_point) ** 2
    
    def displacement_quartiles(self, reference_point=()):
        
        distances = []
        
        if reference_point == (): 
            for i, point1 in enumerate(self.lat_long_coordinates):
                for j, point2 in enumerate(self.lat_long_coordinates):
                    if i != j and i < j:
                        distances.append(self.__haversine(point1, point2))
        else:    
            for lat_long_coordinate in self.lat_long_coordinates:
                distances.append(self.__haversine(reference_point, lat_long_coordinate))
                
        minimum = np.min(distances)
        q1 = np.quantile(distances, 0.25)
        median = np.median(distances)
        q3 = np.quantile(distances, 0.75)
        maximum = np.max(distances)
        iqr = q3 - q1
        dis_range = maximum - minimum
        
        return {"MIN": minimum, "Q1": q1, "MEDIAN": median, "Q3": q3, "MAX": maximum, 
                "IQR": iqr, "RANGE": dis_range}
        
    def geo_midpoint(self):
        x_coords = []
        y_coords = []
        z_coords = []
        
        for lat_long_coordinate in self.lat_long_coordinates:
            cartesian_coord = self.__latlong_to_cartesian(lat_long_coordinate)
            x_coords.append(cartesian_coord[0])
            y_coords.append(cartesian_coord[1])
            z_coords.append(cartesian_coord[2])
        
        return self.__cartesian_to_latlong((np.mean(x_coords), np.mean(y_coords), np.mean(z_coords)))
    
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
        return 6371 * c * 1000
    
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
    

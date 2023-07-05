# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 19:41:04 2023

Title: projection.py
Last Updated: GeoJikuu v0.22.1

Description:
This module contains a class for projecting WGS84 coordinates to linear systems.

Please refer to the official documentation for more information.

Author: Kaine Usher (kaine.usher1@gmail.com)
License: Apache 2.0 (see LICENSE file for details)

Copyright (c) 2023, Kaine Usher.
"""
import math
import random
import pyproj

class CartesianProjector:
    
    def __init__(self, project_from):
        accepted_projections = ["wgs84"]
        
        if project_from.lower() in accepted_projections:
            self.__project_from = project_from.lower()
        else:
            raise ValueError("Construction aborted. The following projection is not available: " + 
                             str(project_from) + " to Cartesian")
            
    def project(self, input_coordinates):
        
        if self.__project_from == "wgs84":
            cartesian_coordinates = []
            unit_conversion = 0
            
            for coordinate_tuple in input_coordinates:
                lat = self.__degrees_to_rads(coordinate_tuple[0])
                lon = self.__degrees_to_rads(coordinate_tuple[1])
                
                x = math.cos(lat) * math.cos(lon)
                y = math.cos(lat) * math.sin(lon)
                z = math.sin(lat)
                
                cartesian_coordinates.append((x, y, z))
            
            unit_conversion = self.__calculate_unit_conversion(input_coordinates, cartesian_coordinates)
            cartesian_coordinates = {"cartesian_coordinates": cartesian_coordinates, "unit_conversion": unit_conversion}
            
        return cartesian_coordinates
    
    def inverse_project(self, input_coordinates):
        wgs84_coordinates = []
        
        for coordinate_tuple in input_coordinates:
            x = coordinate_tuple[0]
            y = coordinate_tuple[1]
            z = coordinate_tuple[2]   
            
            lon = math.atan2(y, x)
            hyp = math.sqrt(x * x + y * y)
            lat = math.atan2(z, hyp)
            
            wgs84_coordinates.append((self.__rads_to_degrees(lat), self.__rads_to_degrees(lon)))
        
        return wgs84_coordinates

            
    def __degrees_to_rads(self, value):
        return (value * math.pi) / 180
    
    def __rads_to_degrees(self, value):
        return (value * 180) / math.pi
            
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
    
    def __euclidean_distance(self, x, y):

        euclid_distance = 0
    
        for i in range(0, len(x)):
            euclid_distance += (float(x[i]) - float(y[i]))**2
    
        return euclid_distance**0.5
    
    # Note: This function only deals with WGS84/Cartesian units. More unit conversions will need to be added later.
    def __calculate_unit_conversion(self, input_coordinates, cartesian_coordinates):
        unit_conversion_samples = []
        
        while len(unit_conversion_samples) < math.ceil((len(input_coordinates)*0.1)):
           random_index_one = random.randint(0, len(input_coordinates)-1)
           random_index_two = random.randint(0, len(input_coordinates)-1)
           if random_index_one == random_index_two:
               continue
           
           metres_diff = self.__haversine(input_coordinates[random_index_one], input_coordinates[random_index_two])
           cartesian_diff = self.__euclidean_distance(cartesian_coordinates[random_index_one], cartesian_coordinates[random_index_two])
           unit_conversion_samples.append(metres_diff / cartesian_diff)
            
        unit_conversion_sum = 0
        
        for sample in unit_conversion_samples:
            unit_conversion_sum += sample
            
        return unit_conversion_sum / len(unit_conversion_samples)
            
class MGA2020Projector:
    
    def __init__(self, project_from):
        accepted_projections = ["wgs84"]
        
        if project_from.lower() in accepted_projections:
            self.__project_from = project_from.lower()
        else:
            raise ValueError("Construction aborted. The following projection is not available: " + 
                             str(project_from) + " to MGA2020")
        
    def project(self, input_coordinates):
        
        if self.__project_from == "wgs84":
            mga2020_coordinates = []
            unit_conversion = 0
            
            wgs84 = pyproj.CRS('EPSG:4326')  # WGS84 (latitude/longitude)
            mga2020 = pyproj.CRS('EPSG:7856')  # MGA2020 (projected coordinate system)
            
            transformer = pyproj.Transformer.from_crs(wgs84, mga2020, always_xy=True)
            
            for coordinate in input_coordinates:
                
                easting, northing = transformer.transform(coordinate[1], coordinate[0])
                mga2020_coordinates.append((northing, easting))

            unit_conversion = self.__calculate_unit_conversion(input_coordinates, mga2020_coordinates)
            mga2020_coordinates = {"mga2020_coordinates": mga2020_coordinates, "unit_conversion": unit_conversion}
            
        return mga2020_coordinates
    
    def inverse_project(self, input_coordinates):
        wgs84_coordinates = []
        
        wgs84 = pyproj.CRS('EPSG:4326')  # WGS84 (latitude/longitude)
        mga2020 = pyproj.CRS('EPSG:7856')  # MGA2020 (projected coordinate system)
        
        transformer = pyproj.Transformer.from_crs(mga2020, wgs84, always_xy=True)
        
        for coordinate in input_coordinates:
            longitude, latitude = transformer.transform(coordinate[1], coordinate[0])
            wgs84_coordinates.append((latitude, longitude))
        
        return wgs84_coordinates
        
    def __degrees_to_rads(self, value):
        return (value * math.pi) / 180
    
    def __rads_to_degrees(self, value):
        return (value * 180) / math.pi
            
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
    
    def __euclidean_distance(self, x, y):
        euclid_distance = 0
    
        for i in range(0, len(x)):
            euclid_distance += (float(x[i]) - float(y[i]))**2
    
        return euclid_distance**0.5
    
    def __calculate_unit_conversion(self, input_coordinates, mga2020_coordinates):
        unit_conversion_samples = []
        while len(unit_conversion_samples) < math.ceil((len(input_coordinates)*0.1)):
           random_index_one = random.randint(0, len(input_coordinates)-1)
           random_index_two = random.randint(0, len(input_coordinates)-1)
           if random_index_one == random_index_two:
               continue
           
           metres_diff = self.__haversine(input_coordinates[random_index_one], input_coordinates[random_index_two])
           mga2020_diff = self.__euclidean_distance(mga2020_coordinates[random_index_one], mga2020_coordinates[random_index_two])
           
           if mga2020_diff == 0:
               if metres_diff == 0:
                   diff_ratio = 1
               else: 
                   continue
           else:
               diff_ratio = metres_diff / mga2020_diff
        
           unit_conversion_samples.append(diff_ratio)
           
        unit_conversion_sum = 0
        
        for sample in unit_conversion_samples:
            unit_conversion_sum += sample
            
        return unit_conversion_sum / len(unit_conversion_samples)
        
    
class MGA1994Projector:
    
    def __init__(self, project_from):
        accepted_projections = ["wgs84"]
        
        if project_from.lower() in accepted_projections:
            self.__project_from = project_from.lower()
        else:
            raise ValueError("Construction aborted. The following projection is not available: " + 
                             str(project_from) + " to MGA1994")
        
    def project(self, input_coordinates):
        
        if self.__project_from == "wgs84":
            mga1994_coordinates = []
            unit_conversion = 0
            
            wgs84 = pyproj.CRS('EPSG:4326')  # WGS84 (latitude/longitude)
            mga1994 = pyproj.CRS('EPSG:28355')  # MGA1994 (projected coordinate system)
            
            transformer = pyproj.Transformer.from_crs(wgs84, mga1994, always_xy=True)
            
            for coordinate in input_coordinates:
                
                easting, northing = transformer.transform(coordinate[1], coordinate[0])
                mga1994_coordinates.append((northing, easting))

            unit_conversion = self.__calculate_unit_conversion(input_coordinates, mga1994_coordinates)
            mga1994_coordinates = {"mga1994_coordinates": mga1994_coordinates, "unit_conversion": unit_conversion}
            
        return mga1994_coordinates
    
    def inverse_project(self, input_coordinates):
        wgs84_coordinates = []
        
        wgs84 = pyproj.CRS('EPSG:4326')  # WGS84 (latitude/longitude)
        mga1994 = pyproj.CRS('EPSG:28355')  # MGA1994 (projected coordinate system)
        
        transformer = pyproj.Transformer.from_crs(mga1994, wgs84, always_xy=True)
        
        for coordinate in input_coordinates:
            longitude, latitude = transformer.transform(coordinate[1], coordinate[0])
            wgs84_coordinates.append((latitude, longitude))
        
        return wgs84_coordinates
        
    def __degrees_to_rads(self, value):
        return (value * math.pi) / 180
    
    def __rads_to_degrees(self, value):
        return (value * 180) / math.pi
            
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
    
    def __euclidean_distance(self, x, y):
        euclid_distance = 0
    
        for i in range(0, len(x)):
            euclid_distance += (float(x[i]) - float(y[i]))**2
    
        return euclid_distance**0.5
    
    def __calculate_unit_conversion(self, input_coordinates, mga2020_coordinates):
        unit_conversion_samples = []
        while len(unit_conversion_samples) < math.ceil((len(input_coordinates)*0.1)):
           random_index_one = random.randint(0, len(input_coordinates)-1)
           random_index_two = random.randint(0, len(input_coordinates)-1)
           if random_index_one == random_index_two:
               continue
           
           metres_diff = self.__haversine(input_coordinates[random_index_one], input_coordinates[random_index_two])
           mga2020_diff = self.__euclidean_distance(mga2020_coordinates[random_index_one], mga2020_coordinates[random_index_two])
           
           if mga2020_diff == 0:
               if metres_diff == 0:
                   diff_ratio = 1
               else: 
                   continue
           else:
               diff_ratio = metres_diff / mga2020_diff
        
           unit_conversion_samples.append(diff_ratio)
           
        unit_conversion_sum = 0
        
        for sample in unit_conversion_samples:
            unit_conversion_sum += sample
            
        return unit_conversion_sum / len(unit_conversion_samples)
        
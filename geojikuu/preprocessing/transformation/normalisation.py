# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 20:56:18 2023
Last Updated: GeoJikuu vx.x.x

Description:
This module contains classes for normalising data, such as scaling and other transformations.

Please refer to the official documentation for more information.

Author: Kaine Usher (kaine.usher1@gmail.com)
License: Apache 2.0 (see LICENSE file for details)

Copyright (c) 2023, Kaine Usher.
"""

class MinMaxScaler:
    
    def __init__(self, data, interval=[0,1]):
        self.__data = data
        self.__min = min(data)
        self.__max = max(data)
        self.__lower_interval = interval[0]
        self.__upper_interval = interval[1]
        
        
    def scale(self, values):
        
        if isinstance(values, list):
            scaled = []
            
            for value in values:
                scaled.append(((value - self.__min) * (self.__upper_interval - self.__lower_interval)) / (self.__max - self.__min) + self.__lower_interval)
            return scaled
            
        else:
            value = values
            return ((value - self.__min) * (self.__upper_interval - self.__lower_interval)) / (self.__max - self.__min) + self.__lower_interval
            
        
    def inverse_scale(self, values):
        
        if isinstance(values, list):
            inverse_scaled = []
            
            for value in values:
                inverse_scaled.append((value - self.__lower_interval) * (self.__max - self.__min) / (self.__upper_interval - self.__lower_interval) + self.__min)
            return inverse_scaled
            
        else:
            value = values
            return (value - self.__lower_interval) * (self.__max - self.__min) / (self.__upper_interval - self.__lower_interval) + self.__min
    
    

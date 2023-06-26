# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 21:44:25 2023

Title: conversion_tools.py
Last Updated: GeoJikuu vx.x.x

Description:
This module contains a collection of conversion tools

Please refer to the official documentation for more information.

Author: Kaine Usher (kaine.usher1@gmail.com)
License: Apache 2.0 (see LICENSE file for details)

Copyright (c) 2023, Kaine Usher.
"""

from datetime import *
import math

class DateConvertor:        

    def __init__(self, date_format_in, date_format_out):
        self.__date_format_in = date_format_in
        self.__date_format_out = date_format_out 

    def date_to_days(self, date):
        
        return self.__temporal_difference(date, datetime.strptime("01/01/0001", "%d/%m/%Y").strftime(self.__date_format_in), self.__date_format_in)
    
    
    def days_to_date(self, days):
        
        return (datetime.strptime("01/01/0001", "%d/%m/%Y").date() + timedelta(days=days)).strftime(self.__date_format_out)
    
    def __temporal_difference(self, time1, time2, date_format):
        
        f_date = datetime.strptime(time2, date_format).date()
        l_date = datetime.strptime(time1, date_format).date()
        delta = l_date - f_date   
        
        return (abs(delta.days))
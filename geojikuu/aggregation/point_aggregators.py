# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 11:30:18 2023

Title: point_aggregators.py
Last Updated: GeoJikuu vx.x.x

Description:
This module contains classes for performing aggregating point data

Please refer to the official documentation for more information.

Author: Kaine Usher (kaine.usher1@gmail.com)
License: Apache 2.0 (see LICENSE file for details)

Copyright (c) 2023, Kaine Usher.
"""

import pandas as pd
import numpy as np

class KNearestNeighbours:
    
    def __init__(self, data, coordinate_label):
        self.__data = data.to_dict()
        self.__coordinate_label = coordinate_label
        

    def aggregate(self, k=1, aggregate_type="sum"):
        
        partitions = self.__partition(k)
        
        points = []
        partition_labels = []
        
        coordinates = self.__data[self.__coordinate_label]
        
        for key, value in coordinates.items():
            partition_index = 0
            points.append(value)
            for partition in partitions:
                for edge in partition:
                    if key in edge:
                        partition_labels.append(partition_index)
                        break
                partition_index = partition_index + 1
            
        aggregate_dict = {
            "partition_labels": partition_labels,
            "points": points,
            }
        
        for key, value in self.__data.items():
            data_entry = []
            if key == self.__coordinate_label:
                continue
            for inner_key, inner_value in value.items():
                data_entry.append(inner_value)
            aggregate_dict[key] = data_entry
            
        
        df = pd.DataFrame.from_dict(aggregate_dict).drop("points", axis=1)
        
        if aggregate_type == "sum":
            df = df.groupby(by=["partition_labels"]).sum()
        elif aggregate_type == "mean":
            df = df.groupby(by=["partition_labels"]).mean()
        else:
            return "Aggregate type not supported."
        
        midpoint_dict = {}
        
        for partition_label in list(set(partition_labels)):
            midpoint_dict[partition_label] = []
            
        for i in range(0,len(partition_labels)):
            midpoint_dict[partition_labels[i]].append(points[i])
        
        for key, value in midpoint_dict.items():
            midpoint_dict[key] = str(self.__midpoint(value))
            
        df2 = pd.DataFrame(midpoint_dict, index=[0]).transpose()
        
        df = df.rename_axis('')
        df2 = df2.rename(columns={df2.columns[0]: 'midpoint'})
        
        return pd.merge(df, df2, left_index=True, right_index=True)


    def __midpoint(self, coord_list):
        
        coord_store = []
        midpoint = []
        for coord in coord_list:
            if type(coord) == str:
                coord_string = coord.strip('(').strip(')').split(", ")
                coord_tuple = tuple([float(i) for i in coord_string])
            else:
                coord_tuple = coord
            coord_store.append(coord_tuple)
            
        for i in range(0,len(coord_store[0])):
            midpoint.append([])
            
        for coord in coord_store:
            for i in range(0,len(midpoint)):
                midpoint[i].append(coord[i])
        
        for i in range(0,len(midpoint)):
            midpoint[i] = np.mean(midpoint[i])
        
        return tuple(midpoint)
    
    def __partition(self, k):
        
        distance_matrix = self.__compute_distance_matrix()
        graph = self.__dict_to_list(self.__compute_graph_dict(distance_matrix, k))
        
        partitioned = []
        
        while graph:
            memory = []
            memory.append(graph[0][0])
            memory.append(graph[0][1])
            partition = [graph.pop(0)]
            
            while memory:
                for edge in graph:
                    if memory[0] in edge:
                        partition.append(edge)
                        memory.append(edge[0])
                        memory.append(edge[1])
                        graph.remove(edge)
                memory.remove(memory[0])
            partitioned.append(partition)
        
        return partitioned
        
    def __compute_graph_dict(self, distance_matrix, k):
        graph_dict = {
            "node": [],
            "neighbour": []
        }
        
        edges_added = [] # Keeps track of edges to prevent duplicates
        for row in distance_matrix:
            if row[0] == "-":
                continue
            
            for i in range(k):
            
                argmin = self.__argmin(row)
                if (str(argmin-1) + str(row[0])) in edges_added:
                    continue
                graph_dict["node"].append(row[0])
                graph_dict["neighbour"].append(argmin-1)
                edges_added.append(str(row[0]) + str(argmin-1))
                
                row[argmin] = -1
                
        return graph_dict
    
    # This function is very ad hoc and should be made redundant in a future update
    def __dict_to_list(self, dictionary):
        
        tuple_list = []
        
        nodes = dictionary["node"]
        neighbours = dictionary["neighbour"]
        
        for i in range(0,len(nodes)):
            tuple_list.append((nodes[i], neighbours[i]))
        
        return tuple_list
        
        
    def __argmin(self, row):
                    
        min_value = max(row)
        min_value_index = np.argmax(row)
            
        for i in range(1, len(row)):
            if row[i] == -1:
                continue
            if i == int(row[0])+1:
                continue
            if row[i] <= min_value:
                min_value = row[i]
                min_value_index = i
        
        return min_value_index
            
    
    def __compute_distance_matrix(self):
        distance_matrix = []
        matrix_row = ['-']
        
        coordinates = self.__data[self.__coordinate_label]
        
        for key in coordinates:
            matrix_row.append(key)
            
        distance_matrix.append(matrix_row)
        
        for key_a, value_a in coordinates.items():
            matrix_row = [key_a]
            for key_b, value_b in coordinates.items():
                if type(value_a) == str or type(value_b) == str:
                    # Convert each string to tuple
                    x_string = value_a.strip('(').strip(')').split(", ")
                    y_string = value_b.strip('(').strip(')').split(", ")
                    x = tuple([float(i) for i in x_string])
                    y = tuple([float(i) for i in y_string])
                else:
                    x = value_a
                    y = value_b
                matrix_row.append(self.__euclidean_distance(x, y))
            distance_matrix.append(matrix_row)
            
        return distance_matrix
            

    def __euclidean_distance(self, x, y):

        euclid_distance = 0
    
        for i in range(0, len(x)):
            euclid_distance += (float(x[i]) - float(y[i]))**2
    
        return euclid_distance**0.5
    
class DistanceBased:
    
    def __init__(self, data, coordinate_label):
        self.__data = data.to_dict()
        self.__coordinate_label = coordinate_label
        
    def aggregate(self, distance, aggregate_type="sum"):
        
        partitions = self.__partition(distance)
        
        points = []
        partition_labels = []
        
        coordinates = self.__data[self.__coordinate_label]
        
        for key, value in coordinates.items():
            partition_index = 0
            points.append(value)
            for partition in partitions:
                for edge in partition:
                    if key in edge:
                        partition_labels.append(partition_index)
                        break
                partition_index = partition_index + 1
            
        aggregate_dict = {
            "partition_labels": partition_labels,
            "points": points,
            }
        
        for key, value in self.__data.items():
            data_entry = []
            if key == self.__coordinate_label:
                continue
            for inner_key, inner_value in value.items():
                data_entry.append(inner_value)
            aggregate_dict[key] = data_entry
            
        
        df = pd.DataFrame.from_dict(aggregate_dict).drop("points", axis=1)
        
        if aggregate_type == "sum":
            df = df.groupby(by=["partition_labels"]).sum()
        elif aggregate_type == "mean":
            df = df.groupby(by=["partition_labels"]).mean()
        else:
            return "Aggregate type not supported."
        
        midpoint_dict = {}
        
        for partition_label in list(set(partition_labels)):
            midpoint_dict[partition_label] = []
            
        for i in range(0,len(partition_labels)):
            midpoint_dict[partition_labels[i]].append(points[i])
        
        for key, value in midpoint_dict.items():
            midpoint_dict[key] = str(self.__midpoint(value))
            
        df2 = pd.DataFrame(midpoint_dict, index=[0]).transpose()
        
        df = df.rename_axis('')
        df2 = df2.rename(columns={df2.columns[0]: 'midpoint'})
        
        return pd.merge(df, df2, left_index=True, right_index=True)
        
    
    def __partition(self, distance):
        
        distance_matrix = self.__compute_distance_matrix()
        graph = self.__dict_to_list(self.__compute_graph_dict(distance_matrix, distance))
        
        partitioned = []
        
        while graph:
            memory = []
            memory.append(graph[0][0])
            memory.append(graph[0][1])
            partition = [graph.pop(0)]
            
            while memory:
                for edge in graph:
                    if memory[0] in edge:
                        partition.append(edge)
                        memory.append(edge[0])
                        memory.append(edge[1])
                        graph.remove(edge)
                memory.remove(memory[0])
            partitioned.append(partition)
        
        return partitioned
    
    def __compute_distance_matrix(self):
        distance_matrix = []
        matrix_row = ['-']
        
        coordinates = self.__data[self.__coordinate_label]
        
        for key in coordinates:
            matrix_row.append(key)
            
        distance_matrix.append(matrix_row)
        
        for key_a, value_a in coordinates.items():
            matrix_row = [key_a]
            for key_b, value_b in coordinates.items():
                if type(value_a) == str or type(value_b) == str:
                    # Convert each string to tuple
                    x_string = value_a.strip('(').strip(')').split(", ")
                    y_string = value_b.strip('(').strip(')').split(", ")
                    x = tuple([float(i) for i in x_string])
                    y = tuple([float(i) for i in y_string])
                else:
                    x = value_a
                    y = value_b
                matrix_row.append(self.__euclidean_distance(x, y))
            distance_matrix.append(matrix_row)
            
        return distance_matrix
    
    def __compute_graph_dict(self, distance_matrix, distance):
        graph_dict = {
            "node": [],
            "neighbour": []
        }
        
        edges_added = [] # Keeps track of edges to prevent duplicates
        index_pos = 1
        for row in distance_matrix:
            if row[0] == "-":
                continue
            
            for i in range(index_pos, len(row)):
                if row[i] <= distance:
                    graph_dict["node"].append(row[0])
                    graph_dict["neighbour"].append(distance_matrix[0][i])
            index_pos += 1
        
        return graph_dict
    
    def __euclidean_distance(self, x, y):

        euclid_distance = 0
    
        for i in range(0, len(x)):
            euclid_distance += (float(x[i]) - float(y[i]))**2
    
        return euclid_distance**0.5
    
    # This function is very ad hoc and should be made redundant in a future update
    def __dict_to_list(self, dictionary):
        
        tuple_list = []
        
        nodes = dictionary["node"]
        neighbours = dictionary["neighbour"]
        
        for i in range(0,len(nodes)):
            tuple_list.append((nodes[i], neighbours[i]))
        
        return tuple_list
    
    def __midpoint(self, coord_list):
        
        coord_store = []
        midpoint = []
        for coord in coord_list:
            if type(coord) == str:
                coord_string = coord.strip('(').strip(')').split(", ")
                coord_tuple = tuple([float(i) for i in coord_string])
            else:
                coord_tuple = coord
            coord_store.append(coord_tuple)
            
        for i in range(0,len(coord_store[0])):
            midpoint.append([])
            
        for coord in coord_store:
            for i in range(0,len(midpoint)):
                midpoint[i].append(coord[i])
        
        for i in range(0,len(midpoint)):
            midpoint[i] = np.mean(midpoint[i])
        
        return tuple(midpoint)
    
class STDistanceBased:
    
    def __init__(self, data, coordinate_label, time_label):
        self.__data = data.to_dict()
        self.__coordinate_label = coordinate_label
        self.__st_coordinate_label = self.__combine_labels(coordinate_label, time_label)
        
    def aggregate(self, spatial_distance, temporal_distance, aggregate_type="sum"):
        
        partitions = self.__partition(spatial_distance, temporal_distance)
        
        points = []
        partition_labels = []
        
        coordinates = self.__data[self.__st_coordinate_label]
        
        for key, value in coordinates.items():
            partition_index = 0
            points.append(value)
            for partition in partitions:
                for edge in partition:
                    if key in edge:
                        partition_labels.append(partition_index)
                        break
                partition_index = partition_index + 1
        
        aggregate_dict = {
            "partition_labels": partition_labels,
            "points": points,
            }
        
        for key, value in self.__data.items():
            data_entry = []
            if key == self.__st_coordinate_label:
                continue
            for inner_key, inner_value in value.items():
                data_entry.append(inner_value)
            aggregate_dict[key] = data_entry
            
        df = pd.DataFrame.from_dict(aggregate_dict).drop("points", axis=1)
        
        if aggregate_type == "sum":
            df = df.groupby(by=["partition_labels"]).sum()
        elif aggregate_type == "mean":
            df = df.groupby(by=["partition_labels"]).mean()
        else:
            return "Aggregate type not supported."
        
        midpoint_dict = {}
        
        for partition_label in list(set(partition_labels)):
            midpoint_dict[partition_label] = []
            
        for i in range(0,len(partition_labels)):
            midpoint_dict[partition_labels[i]].append(points[i])
        
        for key, value in midpoint_dict.items():
            midpoint_dict[key] = str(self.__midpoint(value))
            
        df2 = pd.DataFrame(midpoint_dict, index=[0]).transpose()
        
        df = df.rename_axis('')
        df2 = df2.rename(columns={df2.columns[0]: 'midpoint'})
        
        return pd.merge(df, df2, left_index=True, right_index=True)
    
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
        self.__data["st_temporal_coordinates"] = res_dct
        return "st_temporal_coordinates"

    def __partition(self, spatial_distance, temporal_distance):
        
        distance_matrix = self.__compute_distance_matrix()
        graph = self.__dict_to_list(self.__compute_graph_dict(distance_matrix, spatial_distance, temporal_distance))

        partitioned = []
        
        while graph:
            memory = []
            memory.append(graph[0][0])
            memory.append(graph[0][1])
            partition = [graph.pop(0)]
            
            while memory:
                for edge in graph:
                    if memory[0] in edge:
                        partition.append(edge)
                        memory.append(edge[0])
                        memory.append(edge[1])
                        graph.remove(edge)
                memory.remove(memory[0])
            partitioned.append(partition)
        
        return partitioned
        
    def __compute_distance_matrix(self):
        
        # The distance matrix will need nested tuples (x, y) where x is the spatial distance and y is the
        # temporal difference
        
        distance_matrix = []
        matrix_row = ['-']
        
        coordinates = self.__data[self.__st_coordinate_label]
        
        for key in coordinates:
            matrix_row.append(key)
            
        distance_matrix.append(matrix_row)
        
        for key_a, value_a in coordinates.items():
            matrix_row = [key_a]
            for key_b, value_b in coordinates.items():
                if type(value_a) == str or type(value_b) == str:
                    # Convert each string to tuple
                    x_string = value_a.strip('(').strip(')').split(", ")
                    y_string = value_b.strip('(').strip(')').split(", ")
                    x = tuple([float(i) for i in x_string])
                    y = tuple([float(i) for i in y_string])
                else:
                    x = value_a
                    y = value_b
                matrix_row.append(self.__euclidean_distance(x, y))
            distance_matrix.append(matrix_row)
            
        return distance_matrix
    
    # This is a variation of Euclidean distance that returns a tuple of the (x, y) where x is the spatial distance
    # and y is the temporal difference
    def __euclidean_distance(self, x, y):

        # 1. Work out the number of spatial variables in x and y
        spatial_num = len(self.__data[self.__coordinate_label][0])
        
        # 2. Work out which index represents the temporal variable in x and y
        temporal_index = spatial_num
        
        # 3. Find the Euclidean difference between the spatial variables
        spatial_euclid_distance = 0
    
        for i in range(0, spatial_num):
            spatial_euclid_distance += (float(x[i]) - float(y[i]))**2
            
        spatial_euclid_distance = spatial_euclid_distance**0.5
        
        # 4. Find the Euclidean difference btween the temporal values 
        #    (really just abs(x_temp - y_temp) since temporal variables are one dimensional)
        temporal_euclid_distance = abs(x[temporal_index] - y[temporal_index])
        
        # 5. Return the results in a tuple
        return (spatial_euclid_distance, temporal_euclid_distance)
        
    # This function is very ad hoc and should be made redundant in a future update
    def __dict_to_list(self, dictionary):
        
        tuple_list = []
        
        nodes = dictionary["node"]
        neighbours = dictionary["neighbour"]
        
        for i in range(0,len(nodes)):
            tuple_list.append((nodes[i], neighbours[i]))
        
        return tuple_list

    def __compute_graph_dict(self, distance_matrix, spatial_distance, temporal_distance):
        graph_dict = {
            "node": [],
            "neighbour": []
        }
        
        edges_added = [] # Keeps track of edges to prevent duplicates
        index_pos = 1
        for row in distance_matrix:
            if row[0] == "-":
                continue
            
            for i in range(index_pos, len(row)):
                if row[i][0] <= spatial_distance and row[i][1] <= temporal_distance:
                    graph_dict["node"].append(row[0])
                    graph_dict["neighbour"].append(distance_matrix[0][i])
            index_pos += 1
            
        # Check for loose points and add them to themselves.
                
        return graph_dict
    
    def __midpoint(self, coord_list):
        
        coord_store = []
        midpoint = []
        for coord in coord_list:
            if type(coord) == str:
                coord_string = coord.strip('(').strip(')').split(", ")
                coord_tuple = tuple([float(i) for i in coord_string])
            else:
                coord_tuple = coord
            coord_store.append(coord_tuple)
            
        for i in range(0,len(coord_store[0])):
            midpoint.append([])
            
        for coord in coord_store:
            for i in range(0,len(midpoint)):
                midpoint[i].append(coord[i])
        
        for i in range(0,len(midpoint)):
            midpoint[i] = np.mean(midpoint[i])
        
        return tuple(midpoint)
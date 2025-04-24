import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import MultiPolygon
import re
from pathlib import Path
import json

# Extract the generator information
def Extract_gen_info(CATS_gens_dir):
    # Data for the generators multiple units
    gen_df = pd.read_csv(CATS_gens_dir)  

    # The bus IDs
    all_gen_buses = np.array(gen_df['bus'])
    all_gen_buses_lon = np.array(gen_df['Lon'])
    all_gen_buses_lat = np.array(gen_df['Lat'])

    # Only keep the unique ID numbers
    gen_buses = np.unique(all_gen_buses)

    gen_buses_lon = []
    gen_buses_lat = []
    gen_buses_units = []

    for gb in gen_buses:
        index = np.where(all_gen_buses == gb)[0][0]
        gen_buses_units.append(len(np.where(all_gen_buses == gb)[0]))
        gen_buses_lon.append(all_gen_buses_lon[index])
        gen_buses_lat.append(all_gen_buses_lat[index])

    # Number of units for each generator
    gen_buses_units = np.array(gen_buses_units)
    # The generator longitude coordinates
    gen_buses_lon = np.array(gen_buses_lon)
    # The generator latitude coordinates
    gen_buses_lat = np.array(gen_buses_lat)

    return gen_buses, gen_buses_units, gen_buses_lon, gen_buses_lat

# Extract the lines information 
def Extract_lines_info(CATS_lines_dir):
    # Load JSON file
    with open(CATS_lines_dir, "r") as f:
        cats_lines = json.load(f)

    # Extract features for the lines and subsections of the lines
    features = cats_lines["features"]

    lines_coords = []
    # Iterate through each feature to file coordinates
    for feature in features:
        # Access the geometry dictionary
        geometry = feature["geometry"]    
        # Access the coordinates
        coords = geometry["coordinates"]    
        # If the coordinates are nested (multi-line), flatten them
        if isinstance(coords[0][0], list):
            lines_coords.append([coord for sublist in coords for coord in sublist])
        else:
            # If not nested, append as is
            lines_coords.append(coords)

    # ID coordinates that are an array of array 
    weird_coord_structure = []
    for i, feature in enumerate(features):
        coordinates = feature["geometry"]["coordinates"]
        longitudes = [coord[0] for coord in coordinates]
        latitudes = [coord[1] for coord in coordinates]

        # Check if first entry is a float
        line_lat = latitudes[0]  
        if not isinstance(line_lat, float):  
            weird_coord_structure.append(i)  
    
    return features, lines_coords, weird_coord_structure

# Find the shp files that correspond to the given percentiles
def ID_WF_files(WF_dir, percentiles):
    # folder for the wildfire directory
    folder_path = Path(WF_dir) 
    if not isinstance(percentiles, list): 
        percentiles = [percentiles]

    # Convert x_values to a regex pattern (format as two-digit numbers)
    x_pattern = "|".join(f"{x:02d}" for x in percentiles)
    pattern = re.compile(fr"\w+-\w+_\d+_\d+_\d+_\d+_({x_pattern})_elm\.shp")

    # Find all .shp files and filter based on the pattern
    shapefiles = [file for file in folder_path.glob("*.shp") if pattern.match(file.name)]
    
    return shapefiles

# ID the lines with multiple circuits and the number of circuits
def ID_multicircuits_lines(CATS_lines_circuits_dir):    
    # Data for the lines for the multiple circuits
    line_circuits_df = pd.read_csv(CATS_lines_circuits_dir)  
    multi_circuit_index = np.array(line_circuits_df['Line_index'])
    multi_circuit_Circuits = np.array(line_circuits_df['Circuits'])

    return multi_circuit_index, multi_circuit_Circuits

# ID off buses and their indices, off generators, buses that lost lines and how many were lost
def ID_off_and_lost_info(bus_info, WF_lines_index, gen_buses, display_info = False):
    # ID the nodes that lost all/some of their lines.
    off_nodes_index = []
    off_nodes = []
    off_gen_nodes = []
    lost_lines_index = []
    lost_lines = []

    for i in range(len(bus_info['id'])):
        # The lines not affected by the wildfire that are connected to bus i
        not_WF_lines = np.setdiff1d(bus_info['bus_lines_index'][i], WF_lines_index)
        # If all the lines connected to bus i are affected by the wildfire
        if len(not_WF_lines) == 0:
            off_nodes_index.append(i)
            off_nodes.append(bus_info['id'][i])
            # Check if  bus is a generator of load
            if bus_info['id'][i] in gen_buses:
                off_gen_nodes.append(bus_info['id'][i])
                if display_info:
                    print('Bus (generator) ' + str(bus_info['id'][i]) + ' is off.')
            else:
                if display_info:
                    print('Bus (load) ' + str(bus_info['id'][i]) + ' is off.')
        # The number of lines not affected by the wildfire are less than the total number of lines
        # connected to bus i
        elif len(not_WF_lines) < bus_info['total_lines'][i]:
            # The number of lines not affected by the wildfire
            lines_left = bus_info['total_lines'][i] - len(not_WF_lines)
            lost_lines_index.append(i)
            lost_lines.append(lines_left)
            if bus_info['id'][i] in gen_buses:
                if display_info:
                    print('Bus (generator) ' + str(bus_info['id'][i]) + ' lost ' + str(lines_left) + ' out of ' + str(bus_info['total_lines'][i]))
            else:
                if display_info:
                    print('Bus (load) ' + str(bus_info['id'][i]) + ' lost ' + str(lines_left) + ' out of ' + str(bus_info['total_lines'][i]))

    return np.array(off_nodes_index), np.array(off_nodes), np.array(off_gen_nodes), np.array(lost_lines_index), np.array(lost_lines)

# ID bus information
def ID_bus_info(CATS_buses_dir, features):
    bus_data = pd.read_csv(CATS_buses_dir) 

    # The from and to bus for each line
    from_bus, to_bus = endpoint_info(features)

    # extract the longitude and latitude coordinates
    nodes_lat = bus_data.iloc[:, 4].values
    nodes_lon = bus_data.iloc[:, 5].values

    bus_id = np.array(bus_data['bus_i'])

    bus_lines_index = []
    total_bus_lines_index = []
    for x in bus_id:
        # Get indices where the bus id appears in the set of from and to bus
        indices1 = np.where(from_bus == x)[0]  
        indices2 = np.where(to_bus == x)[0]  

        # Combine indices into a single 1D array
        combined_indices = np.concatenate((indices1, indices2))
        bus_lines_index.append(combined_indices)
        total_bus_lines_index.append(len(combined_indices))

    # Dictornary containg the bus id, the index of lines connected to the bus and total number of line
    bus_info = {'id':bus_id, 'bus_lines_index': bus_lines_index, 'total_lines': np.array(total_bus_lines_index)}

    return bus_info, nodes_lat, nodes_lon

# ID the transformers indices
def ID_transformer(CATS_transformer_dir):    
    # Data for the transformer
    transformer_df = pd.read_csv(CATS_transformer_dir)  
    transformer_index = np.array(transformer_df['Transformer_index'])

    return transformer_index

# Extract the wildfire lontitude and latitude coordinates
def Extract_WF_Coordinates(fire_area_path):

    # Extract wildfire coordinates
    fire_area_gdf = gpd.read_file(fire_area_path)
    # Map coordinates to correspond to California
    fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326) #EPSG 4326 is for California

    # Extract all points from all MultiPolygons
    all_points = []
    for idx, row in fire_area_4326_gdf.iterrows():
        geom = row["geometry"]
        if isinstance(geom, MultiPolygon):
            for polygon in geom.geoms:
                all_points.extend(list(polygon.exterior.coords))
        elif isinstance(geom, Polygon):  # In case there are single polygons
            all_points.extend(list(geom.exterior.coords))

    # Convert points to a GeoDataFrame
    points_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(*zip(*all_points)), crs=fire_area_gdf.crs)

    # Extract longitude and latitude separately
    WF_lon = np.array(points_gdf.geometry.x)
    WF_lat = np.array(points_gdf.geometry.y)

    return WF_lon, WF_lat

# Find the generator within the safe distance and and the number of units correspond to the generator
# that will be added to the contingency file
def ID_Gen_near_WF(gen_buses_units, gen_buses, gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, SD):

    # ID the generators within the safe distance
    out_gen_indices, out_gen_indices_dis = generators_within_safe_distance(gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, SD)
    out_gen_units = gen_buses_units[out_gen_indices]

    contin_gen_unit = []
    contin_gen_bus = []
    for i in range(len(out_gen_units)):
        for k in range(out_gen_units[i]):
            contin_gen_unit.append(k+1)
            contin_gen_bus.append(gen_buses[out_gen_indices[i]])
    # The generator and the number of units that will go into the contingency file
    contin_gen_unit = np.array(contin_gen_unit)
    contin_gen_bus = np.array(contin_gen_bus)

    return contin_gen_unit, contin_gen_bus

# Same as ID_Gen_near_WF but now the indices of the out generators are given
def ID_Gen_Contin_info(gen_buses, gen_buses_units, out_gen_indices):
    
    out_gen_units = gen_buses_units[out_gen_indices]

    contin_gen_unit = []
    contin_gen_bus = []
    for i in range(len(out_gen_units)):
        for k in range(out_gen_units[i]):
            contin_gen_unit.append(k+1)
            contin_gen_bus.append(gen_buses[out_gen_indices[i]])
    contin_gen_unit = np.array(contin_gen_unit)
    contin_gen_bus = np.array(contin_gen_bus)

    return contin_gen_unit, contin_gen_bus

# ID the lines withing the safe distance that will be added to the contingency file
def ID_Lines_near_WF(features, WF_lon, WF_lat, SD, weird_coord_structure, multi_circuit_index, multi_circuit_Circuits):
    # ID the lines with the safe distance
    WF_lines_index, WF_lines_id, WF_line_dis = lines_within_safe_distance(features, WF_lon, WF_lat, SD, weird_coord_structure)

    contin_lines_from_buses = []
    contin_lines_to_buses = []
    contin_circuit = []

    for i in WF_lines_index:
        if i in multi_circuit_index:
            index = np.where(multi_circuit_index == i)[0][0]
            for k in range(multi_circuit_Circuits[index]):
                contin_lines_from_buses.append(features[i]['properties']['f_bus'])
                contin_lines_to_buses.append(features[i]['properties']['t_bus'])
                contin_circuit.append(k + 1)
        else:
                contin_lines_from_buses.append(features[i]['properties']['f_bus'])
                contin_lines_to_buses.append(features[i]['properties']['t_bus'])
                contin_circuit.append(1)

    # The lines within the safe distance from, to bus and the number of circuits
    contin_lines_from_buses = np.array(contin_lines_from_buses)
    contin_lines_to_buses = np.array(contin_lines_to_buses)
    contin_circuit = np.array(contin_circuit)

    return contin_lines_from_buses, contin_lines_to_buses, contin_circuit

# Same as ID_Lines_near_WF but the indices of the lines within the safe distance are given
def ID_Lines_Contin_info(features, WF_lines_index, multi_circuit_index, multi_circuit_Circuits):
    # ID the lines with the safe distance
    contin_lines_from_buses = []
    contin_lines_to_buses = []
    contin_circuit = []

    for i in WF_lines_index:
        if i in multi_circuit_index:
            index = np.where(multi_circuit_index == i)[0][0]
            for k in range(multi_circuit_Circuits[index]):
                contin_lines_from_buses.append(features[i]['properties']['f_bus'])
                contin_lines_to_buses.append(features[i]['properties']['t_bus'])
                contin_circuit.append(k + 1)
        else:
                contin_lines_from_buses.append(features[i]['properties']['f_bus'])
                contin_lines_to_buses.append(features[i]['properties']['t_bus'])
                contin_circuit.append(1)

    contin_lines_from_buses = np.array(contin_lines_from_buses)
    contin_lines_to_buses = np.array(contin_lines_to_buses)
    contin_circuit = np.array(contin_circuit)

    return contin_lines_from_buses, contin_lines_to_buses, contin_circuit

# Function to calculate distance between two points using the Haversine formula
def haversine_distance(lat1, lon1, lat2, lon2):
    # Earth radius in miles
    r = 3959.0

    # Convert degrees to radians
    deg2rad = np.pi / 180
    lat1  = lat1 * deg2rad 
    lon1  = lon1 * deg2rad
    lat2  = lat2 * deg2rad
    lon2  = lon2 * deg2rad

    # Differences in coordinates
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Distance in miles
    return r * c  # Return distance in miles

# Find from and to bus for each line
def endpoint_info(features):
    from_bus = []
    to_bus = []
    for i in range(len(features)):
        from_bus.append(features[i]['properties']['f_bus'])
        to_bus.append(features[i]['properties']['t_bus'])
    return np.array(from_bus), np.array(to_bus)

# Vectorized version of the Haversine function for array inputs
def haversine_distance_vectorized(lat1, lon1, lat2, lon2):
    # Earth radius in miles
    r = 3959.0

    lat1 = np.array(lat1)
    lon1 = np.array(lon1)
    lat2 = np.array(lat2)
    lon2 = np.array(lon2)

    # Convert degrees to radians
    deg2rad = np.pi / 180
    lat1  = lat1 * deg2rad 
    lon1  = lon1 * deg2rad
    lat2  = lat2 * deg2rad
    lon2  = lon2 * deg2rad

    # Differences in coordinates
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    # c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    # Distance in miles
    return r * c

# Extract the longitude and latitude coordinate from the set of lines that are defined as a set of set instead
# of a set of points
def ID_lines_coords(features, weird_coord_structure):
    lon = []
    lat = []

    for i in range(len(features)):
        if i in weird_coord_structure:  
            longitudes = []
            latitudes = []
            for l in range(len(features[i]["geometry"]["coordinates"])):
                longitudes.extend([coord[0] for coord in features[i]["geometry"]["coordinates"][l]])
                latitudes.extend([coord[1] for coord in features[i]["geometry"]["coordinates"][l]])

        else:
            longitudes = [coord[0] for coord in features[i]["geometry"]["coordinates"]]
            latitudes = [coord[1] for coord in features[i]["geometry"]["coordinates"]]

            lon.append(longitudes)
            lat.append(latitudes)

    return lon, lat

# Calculate the nearest point on the line to the fire perimeter and ID if they are within SD miles
def lines_within_safe_distance(features, WF_lon, WF_lat, SD, weird_coord_structure):
    indices = []
    ids = []
    indices_dis = []

    for i in range(len(features)):
        dis = []
        if i in weird_coord_structure:  
            for l in range(len(features[i]["geometry"]["coordinates"])):
                longitudes = [coord[0] for coord in features[i]["geometry"]["coordinates"][l]]
                latitudes  = [coord[1] for coord in features[i]["geometry"]["coordinates"][l]]

                for k in range(len(longitudes)):
                    line_lat = latitudes[k]
                    line_lon = longitudes[k]
                    dis.append(haversine_distance_vectorized(WF_lat, WF_lon, line_lat, line_lon))

        else:
            longitudes = [coord[0] for coord in features[i]["geometry"]["coordinates"]]
            latitudes = [coord[1] for coord in features[i]["geometry"]["coordinates"]]
            for k in range(len(longitudes)):
                line_lat = latitudes[k]
                line_lon = longitudes[k]
                dis.append(haversine_distance_vectorized(WF_lat, WF_lon, line_lat, line_lon))

        min_dis = np.min(dis)
        if min_dis < SD:
            indices_dis.append(min_dis)
            indices.append(i)
            ids.append(features[i]["id"])

    return np.array(indices), np.array(ids), np.array(indices_dis)

# Calculate the distance between the generators and wildfire perimeter then ID if they are within SD miles.
def generators_within_safe_distance(gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, SD):
    indices = []
    indices_dis = []

    for i in range(len(gen_buses_lon)):
        dis = np.min(haversine_distance_vectorized(WF_lat, WF_lon, gen_buses_lat[i], gen_buses_lon[i]))

        if dis < SD:
            indices_dis.append(dis)
            indices.append(i)

    return np.array(indices), np.array(indices_dis)
import os
from required_funtions import *

# Indicates what type of wildfire data (historic or probabilistic) 
# If historic (json), then output the different titles
# If probabilisic (shp), then output the different percentile
def extract_WF_dir_info(directory, print_con_file_details = True):
    json_files = []
    shp_files = []

    id_type = False

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.json'):
                if id_type == False:
                    file_type = 'json'
                    id_type = True

                json_files.append(os.path.join(root, file))
            elif file.endswith('.shp'):
                if id_type == False:
                    file_type = 'shp'
                    id_type = True

                shp_files.append(os.path.join(root, file))

    if file_type == 'json':
        # Extract the json file names
        filenames = [os.path.basename(f) for f in json_files]
        print("The different json files:")
        print(filenames)
        print("\n")
        print("To create the con file use function: \n create_WF_con_file_json")
        return file_type, filenames

    else:
        # Extract the percentile from the files
        pattern = re.compile(r'_(\d+)_elm\.shp$')
        percentiles = np.sort([int(pattern.search(f).group(1)) for f in shp_files if pattern.search(f)])
        print("The different percentile from the shp file:")
        print(percentiles)
        if print_con_file_details:
            print("\n")
            print("To create the con file use function: \n create_WF_con_file_shp")
        return file_type, percentiles

# Check if the set of percentile chosen are in the set of possible percentile
# The possible percentile come from the function extract_WF_dir_info
def ID_percentiles_in_WF_dir(percentiles, possible_percentiles):
    if len(percentiles) != len(np.intersect1d(percentiles, possible_percentiles)):
        print('\nThe following percentile/s: ' + str(np.setdiff1d(percentiles, possible_percentiles)) + ' are not one of the possible percentile')
        percentiles = np.intersect1d(percentiles, possible_percentiles)
        return list(percentiles)
    else:
        print("\nAll the percentiles belong in the set of possible percentiles")
        return percentiles

# Extract the wildfire coordinates of the single json file
def extract_WF_coords_json(WF_dir, json_filename):
    WF_path = WF_dir + '/' + json_filename
    WF_data_gdf = gpd.read_file(WF_path)
    # Extract all points from all MultiPolygons
    all_points = []
    for idx, row in WF_data_gdf.iterrows():
        geom = row["geometry"]
        if isinstance(geom, MultiPolygon):
            for polygon in geom.geoms:
                all_points.extend(list(polygon.exterior.coords))
        elif isinstance(geom, Polygon):  # In case there are single polygons
            all_points.extend(list(geom.exterior.coords))

    # Convert points to a GeoDataFrame
    points_gdf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(*zip(*all_points)), crs=WF_data_gdf.crs)

    # Extract longitude and latitude separately
    WF_lon = np.array(points_gdf.geometry.x)
    WF_lat = np.array(points_gdf.geometry.y)

    return WF_data_gdf, WF_lon, WF_lat

# Extract the wildfire coordinates of the single percentile
def extract_WF_coords_shp(WF_dir, possible_percentiles, percentile):
    if percentile not in possible_percentiles:
        raise ValueError(
            str(percentile) + " is not a possible percentile. Valid options are: " + ", ".join(map(str, possible_percentiles))
        )
    shapefiles = ID_WF_files(WF_dir, percentile)
    fire_area_path = WF_dir + "/" + shapefiles[0].name
    fire_area_gdf = gpd.read_file(fire_area_path)
    WF_lon, WF_lat = Extract_WF_Coordinates(fire_area_path)

    return fire_area_gdf, WF_lon, WF_lat

# ID the information of element turned off
# Lines and transformer: their ids and the number of circuits
# Generators: their index and the number of units
def ID_WF_Contingency_info_json(WF_dir, CATS_lines_dir, CATS_lines_circuits_dir, CATS_gens_dir, CATS_transformer_dir, SD, json_filenames):

    transformer_index = ID_transformer(CATS_transformer_dir)

    multi_circuit_index, multi_circuit_Circuits = ID_multicircuits_lines(CATS_lines_circuits_dir)

    gen_buses, gen_buses_units, gen_buses_lon, gen_buses_lat = Extract_gen_info(CATS_gens_dir)

    features, lines_coords, weird_coord_structure = Extract_lines_info(CATS_lines_dir)

    # Information indicating number of lines, transformers, and generators turned off
    off_info = []

    # Information needed for the contingency file
    contin_unit = []
    contin_gen = []
    contin_from = []
    contin_to = []
    contin_cir = []

    # Information needed for the plots
    WF_data_gdfs = []
    WF_lons = []
    WF_lats = []
    out_gen_indixs = []
    WF_lines_indexs = []

    # iterate through the json files and extract the desired information
    for wf in json_filenames:

        WF_data_gdf, WF_lon, WF_lat = extract_WF_coords_json(WF_dir, wf)

        out_gen_indix, out_gen_indices_dis = generators_within_safe_distance(gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, SD)
        
        if len(out_gen_indices_dis) > 0:
            contin_gen_unit, contin_gen_bus = ID_Gen_Contin_info(gen_buses, gen_buses_units, out_gen_indix)
        
        WF_lines_index, WF_lines_id, WF_line_dis = lines_within_safe_distance(features, WF_lon, WF_lat, SD, weird_coord_structure)
        if len(WF_lines_index) > 0:
            contin_lines_from_buses, contin_lines_to_buses, contin_circuit = ID_Lines_Contin_info(features, WF_lines_index, multi_circuit_index, multi_circuit_Circuits)
        
        WF_trans_index = np.intersect1d(WF_lines_index, transformer_index)
        WF_non_trans_index = np.setdiff1d(WF_lines_index, WF_trans_index)

        off_info.append([len(WF_non_trans_index), len(WF_trans_index), len(out_gen_indix)])

        WF_data_gdfs.append(WF_data_gdf)
        WF_lons.append(WF_lon)
        WF_lats.append(WF_lat)
        out_gen_indixs.append(out_gen_indix)
        WF_lines_indexs.append(WF_lines_index)

        if (len(WF_lines_index) > 0) and (len(out_gen_indix) > 0):
            contin_unit.append(contin_gen_unit)
            contin_gen.append(contin_gen_bus)
            contin_from.append(contin_lines_from_buses)
            contin_to.append(contin_lines_to_buses)
            contin_cir.append(contin_circuit)
        elif len(WF_lines_index) > 0:
            contin_unit.append([])
            contin_gen.append([])
            contin_from.append(contin_lines_from_buses)
            contin_to.append(contin_lines_to_buses)
            contin_cir.append(contin_circuit)
        elif len(out_gen_indix) > 0:
            contin_unit.append(contin_gen_unit)
            contin_gen.append(contin_gen_bus)
            contin_from.append([])
            contin_to.append([])
            contin_cir.append([])
        
    return contin_unit, contin_gen, contin_from, contin_to, contin_cir, WF_data_gdfs, WF_lons, WF_lats, out_gen_indixs, WF_lines_indexs, off_info

# ID the information of element turned off
# Lines and transformer: their ids and the number of circuits
# Generators: their index and the number of units
def ID_WF_Contingency_info_shp(WF_dir, CATS_lines_dir, CATS_lines_circuits_dir, CATS_gens_dir, CATS_transformer_dir, SD, percentiles, possible_percentiles):

    transformer_index = ID_transformer(CATS_transformer_dir)

    multi_circuit_index, multi_circuit_Circuits = ID_multicircuits_lines(CATS_lines_circuits_dir)

    gen_buses, gen_buses_units, gen_buses_lon, gen_buses_lat = Extract_gen_info(CATS_gens_dir)

    features, lines_coords, weird_coord_structure = Extract_lines_info(CATS_lines_dir)

    if len(percentiles) != len(np.intersect1d(percentiles, possible_percentiles)):
        print('The percentile/s ' + str(np.setdiff1d(percentiles, possible_percentiles)) + ' is/are not one of the possible percentile')
        percentiles = np.intersect1d(percentiles, possible_percentiles)
        print('Will create the contingency file for the following percentiles: ' + str(percentiles) + "\n")

    # Information indicating number of lines, transformers, and generators turned off
    off_info = []

    # Information needed for the contingency file
    contin_unit = []
    contin_gen = []
    contin_from = []
    contin_to = []
    contin_cir = []

    # Information needed for the plots
    WF_data_gdfs = []
    WF_lons = []
    WF_lats = []
    out_gen_indixs = []
    WF_lines_indexs = []

    # iterate through the shp files and extract the desired information
    for wf in percentiles:
        WF_data_gdf, WF_lon, WF_lat = extract_WF_coords_shp(WF_dir, possible_percentiles, wf)

        out_gen_indix, out_gen_indices_dis = generators_within_safe_distance(gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, SD)
        if len(out_gen_indices_dis) > 0:
            contin_gen_unit, contin_gen_bus = ID_Gen_Contin_info(gen_buses, gen_buses_units, out_gen_indix)
        
        WF_lines_index, WF_lines_id, WF_line_dis = lines_within_safe_distance(features, WF_lon, WF_lat, SD, weird_coord_structure)
        if len(WF_lines_index) > 0:
            contin_lines_from_buses, contin_lines_to_buses, contin_circuit = ID_Lines_Contin_info(features, WF_lines_index, multi_circuit_index, multi_circuit_Circuits)
        
        WF_trans_index = np.intersect1d(WF_lines_index, transformer_index)
        WF_non_trans_index = np.setdiff1d(WF_lines_index, WF_trans_index)

        off_info.append([len(WF_non_trans_index), len(WF_trans_index), len(out_gen_indix)])

        WF_data_gdfs.append(WF_data_gdf)
        WF_lons.append(WF_lon)
        WF_lats.append(WF_lat)
        out_gen_indixs.append(out_gen_indix)
        WF_lines_indexs.append(WF_lines_index)
        
        if (len(WF_lines_index) > 0) and (len(out_gen_indix) > 0):
            contin_unit.append(contin_gen_unit)
            contin_gen.append(contin_gen_bus)
            contin_from.append(contin_lines_from_buses)
            contin_to.append(contin_lines_to_buses)
            contin_cir.append(contin_circuit)
        elif len(WF_lines_index) > 0:
            contin_unit.append([])
            contin_gen.append([])
            contin_from.append(contin_lines_from_buses)
            contin_to.append(contin_lines_to_buses)
            contin_cir.append(contin_circuit)
        elif len(out_gen_indix) > 0:
            contin_unit.append(contin_gen_unit)
            contin_gen.append(contin_gen_bus)
            contin_from.append([])
            contin_to.append([])
            contin_cir.append([])
        
    return contin_unit, contin_gen, contin_from, contin_to, contin_cir, WF_data_gdfs, WF_lons, WF_lats, out_gen_indixs, WF_lines_indexs, off_info

# Create con file for the different contingency provide.
def create_WF_con_file_json(WF_dir, CATS_lines_dir, CATS_lines_circuits_dir, CATS_gens_dir, CATS_transformer_dir, output_dir, SD, WF_file_type, json_files, con_filename = None):
    # Inform the person if wrong code was used
    if WF_file_type == "shp":
        raise ValueError(
            "The type of data given was .shp instead of .json. Use the following function for this data: create_WF_con_file_shp " 
        )
    
    # extract different json files
    json_filenames = [name[:-5] for name in json_files if name.endswith(".json")]
        
    # All the off elements
    contin_unit, contin_gen, contin_from, contin_to, contin_cir, WF_data_gdfs, WF_lons, WF_lats, out_gen_indixs, WF_lines_indexs, off_info = ID_WF_Contingency_info_json(WF_dir, CATS_lines_dir, CATS_lines_circuits_dir, CATS_gens_dir, CATS_transformer_dir, SD, json_files)

    # Inform the person what cases are finish and the elements turned off
    print("The off information due to the wildfire by congiven contingency:")
    for k in range(len(json_files)):
        print("Contingency: " + json_filenames[k] + " with SD = " + str(SD))
        print("Number of off lines: " + str(off_info[k][0]))
        print("Number of off transformers: " + str(off_info[k][1]))
        print("Number of off generator buses: " + str(off_info[k][2]) + "\n")

    # The process to create the con file
    con_names = []
    if con_filename == None:
        filename = "/case.con"
    else:
        if con_filename[0] != '/':
            filename = '/' + con_filename

    with open(output_dir + filename, "w") as f:
        for k in range(len(json_files)):

            con_names.append(json_filenames[k] + "_SD_" + str(SD))

            f.write("CONTINGENCY "+ json_filenames[k] + "_SD_" + str(SD) + "\n")

            for g in range(len(contin_unit[k])):
                f.write("REMOVE UNIT " + str(contin_unit[k][g]) +" FROM BUS " + str(contin_gen[k][g]) + "\n")

            for b in range(len(contin_cir[k])):
                f.write("OPEN BRANCH FROM BUS " + str(contin_from[k][b] )+" TO BUS " + str(contin_to[k][b])  + " CIRCUIT " + str(contin_cir[k][b]) + "\n")
            
            f.write("END\n")

        f.write("END\n")

    print("Finished creating " + filename)

    return WF_data_gdfs, WF_lons, WF_lats, out_gen_indixs, WF_lines_indexs, con_names

# Create con file for the different contingency provide.
def create_WF_con_file_shp(WF_dir, CATS_lines_dir, CATS_lines_circuits_dir, CATS_gens_dir, CATS_transformer_dir, output_dir, SD, percentiles, WF_file_type, possible_percentiles, con_filename = None):
    # Inform the person if wrong code was used
    if WF_file_type == "json":
        raise ValueError(
            "The type of data given was .json instead of .shp. Use the following function for this data: create_WF_con_file_json " 
        )
        
    # ID what percentiles are possible
    if len(percentiles) != len(np.intersect1d(percentiles, possible_percentiles)):
        print('The following percentile/s: ' + str(np.setdiff1d(percentiles, possible_percentiles)) + ' are not one of the possible percentile')
        percentiles = np.intersect1d(percentiles, possible_percentiles)
        print('Will create the contingency file for the following percentiles: ' + str(percentiles) + "\n")

    # All the off elements
    contin_unit, contin_gen, contin_from, contin_to, contin_cir, WF_data_gdfs, WF_lons, WF_lats, out_gen_indixs, WF_lines_indexs, off_info = ID_WF_Contingency_info_shp(WF_dir, CATS_lines_dir, CATS_lines_circuits_dir, CATS_gens_dir, CATS_transformer_dir, SD, percentiles, possible_percentiles)

    # Inform the person what cases are finish and the elements turned off
    print("The off information due to the wildfire by given contingency:")
    for k in range(len(percentiles)):
        print("Contingency: Percentile_" + str(percentiles[k]) + " with SD = " + str(SD))
        print("Number of off lines: " + str(off_info[k][0]))
        print("Number of off transformers: " + str(off_info[k][1]))
        print("Number of off generator buses: " + str(off_info[k][2]) + "\n")

    # The process to create the con file
    con_names = []
    if con_filename == None:
        filename = "/case.con"
    else:
        if con_filename[0] != '/':
            filename = '/' + con_filename

    with open(output_dir + filename, "w") as f:
        for k in range(len(percentiles)):

            con_names.append("Percentile_"+ str(percentiles[k]) + "_SD_" + str(SD) )

            f.write("CONTINGENCY Percentile_"+ str(percentiles[k]) + "_SD_" + str(SD) + "\n")

            for g in range(len(contin_unit[k])):
                f.write("REMOVE UNIT " + str(contin_unit[k][g]) +" FROM BUS " + str(contin_gen[k][g]) + "\n")

            for b in range(len(contin_cir[k])):
                f.write("OPEN BRANCH FROM BUS " + str(contin_from[k][b] )+" TO BUS " + str(contin_to[k][b])  + " CIRCUIT " + str(contin_cir[k][b]) + "\n")
            
            f.write("END\n")

        f.write("END\n")

    print("Finished creating " + filename)

    return WF_data_gdfs, WF_lons, WF_lats, out_gen_indixs, WF_lines_indexs, con_names

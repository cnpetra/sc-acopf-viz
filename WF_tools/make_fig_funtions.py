from matplotlib.lines import Line2D
import contextily as ctx
import ast  
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re
import os

import matplotlib.lines as mlines
import matplotlib.colors as mcolors

from required_funtions import *
from validation_functions import *

# Make figures that overlay the different wildfire area on top of one another
def Make_WF_fig_Areas(percentiles, WF_dir):
    # Define a colormap and create transparent colors
    colors = plt.cm.plasma(np.linspace(0, 1, len(percentiles)))
    colors_with_alpha = [(c[0], c[1], c[2], 0.5) for c in colors] 

    # The different wildfire by percentiles
    shapefiles = ID_WF_files(WF_dir, percentiles)

    fig, ax = plt.subplots(figsize=(10, 10))

    # Iterate in reverse order so lower percentiles are on top
    legend_handles = []
    for i in reversed(range(len(percentiles))):
        fire_area_path = WF_dir + "/" + shapefiles[i].name
        fire_area_gdf = gpd.read_file(fire_area_path)
        fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326)  # Convert to EPSG 4326

        # Plot each fire area with transparency
        fire_area_4326_gdf.plot(ax=ax, edgecolor='k', color=colors[i], alpha=0.5)

        # Add a custom legend entry
        legend_handles.append(Line2D([0], [0], marker='o', color='w', markersize=10, 
                                    markerfacecolor=colors[i], label=f"Percentile {percentiles[i]}"))

    # Add a basemap
    ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string())

    # Add a legend with percentile labels
    ax.legend(handles=legend_handles, title="Fire Percentiles", loc='upper right')

    # Set title
    plt.title("Wildfire Area by Percentile")

    # Show the plot
    plt.show()

# Make figures that overlay the different wildfire perimeter on top of one another
def Make_WF_fig_Perimeter(percentiles, WF_dir):
   # Define colormap
    colors = plt.cm.plasma(np.linspace(0, 1, len(percentiles)))

    shapefiles = ID_WF_files(WF_dir, percentiles)

    fig, ax = plt.subplots(figsize=(10, 10))
    legend_handles = []

    # Iterate normally to layer properly (newer fires on top)
    for i in reversed(range(len(percentiles))):
        fire_area_path = WF_dir + "/" + shapefiles[i].name
        fire_area_gdf = gpd.read_file(fire_area_path)
        fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326)  # Convert to EPSG 4326

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

        # Plot points with GeoDataFrame's plot function, color based on the current percentile
        points_gdf.plot(ax=ax, color=colors[i], markersize=1, label=f"Percentile {percentiles[i]}")

        # Create a custom legend handle
        legend_handles.append(Line2D([0], [0], marker='o', color='w', markersize=10, 
                                    markerfacecolor=colors[i], markeredgewidth=1, label=f"Percentile {percentiles[i]}"))

    # Add a legend
    ax.legend(handles=legend_handles, title="Fire Percentiles", loc='upper right')

    # Set title
    plt.title("Wildfire Perimeter by Percentile")

    ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string(), zorder=0)

    # Show the plot
    plt.show() 

# Function required for Off_info_vai_SD and Off_info_vai_SD_Percentile
def parse_tuple(s):
    if isinstance(s, str) and re.match(r"^\(\d+,\s*\d+,\s*\d+\)$", s.strip()):
        return ast.literal_eval(s.strip())  
    return (0, 0, 0)  

# ID the number of line/transfomers/generators are turned of when SD is varied for historic data
def Off_info_vai_SD(WF_dir, CATS_lines_dir, CATS_transformer_dir, CATS_gens_dir, SD, save_csv = False, csv_path = None, display_off_info = False):

    transformer_index = ID_transformer(CATS_transformer_dir)

    gen_buses, gen_buses_units, gen_buses_lon, gen_buses_lat = Extract_gen_info(CATS_gens_dir)

    features, lines_coords, weird_coord_structure = Extract_lines_info(CATS_lines_dir)

    df = pd.DataFrame()
    fire_area_gdf = gpd.read_file(WF_dir)
    WF_lon, WF_lat = Extract_WF_Coordinates(WF_dir)

    for i in range(len(SD)):
        sd = SD[i]

        # Off info for generators
        out_gen_indices, out_gen_indices_dis = generators_within_safe_distance(gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, sd)

        # Off info for lines
        WF_lines_index, WF_lines_id, WF_line_dis = lines_within_safe_distance(features, WF_lon, WF_lat, sd, weird_coord_structure)
        
        # Find off lines and off transformers
        WF_trans_index = np.intersect1d(WF_lines_index, transformer_index)
        WF_non_trans_index = np.setdiff1d(WF_lines_index, WF_trans_index)

        print("Finished calculating infomation for  SD: " + str(SD[i]))
        if display_off_info:
            print("The off information due to the wildfire:")   
            print("Number of off lines: " + str(len(WF_non_trans_index)))
            print("Number of off transformers: " + str(len(WF_trans_index)))
            print("Number of off generator buses: " + str(len(out_gen_indices)) + "\n")    

        # Store the data
        data = {
                    str(sd): (len(WF_non_trans_index), len(WF_trans_index), len(out_gen_indices))
                }
    
    # Convert the tuples to strings
    data = {key: [str(val) for val in value] for key, value in data.items()}
    
    if save_csv == True:
        if csv_path == None:
            df.to_csv("./off_info_data.csv", index=True, index_label="Distance / Percentiles")
        else:
            df.to_csv(csv_path, index=True, index_label="Distance / Percentiles")

    return df

# ID the number of line/transfomers/generators are turned of when SD and percentiles are varied for probabilistic data
def Off_info_vai_SD_Percentile(WF_dir, CATS_lines_dir, CATS_transformer_dir, CATS_gens_dir, SD, percentiles, save_csv = False, csv_path = None, display_off_info = False):

    transformer_index = ID_transformer(CATS_transformer_dir)

    gen_buses, gen_buses_units, gen_buses_lon, gen_buses_lat = Extract_gen_info(CATS_gens_dir)

    features, lines_coords, weird_coord_structure = Extract_lines_info(CATS_lines_dir)

    shapefiles = ID_WF_files(WF_dir, percentiles)

    if len(shapefiles) != len(percentiles):
        print("Warning: One of the percentiles data was not within the directory")

    df = pd.DataFrame()
    for wf in range(len(shapefiles)):
        fire_area_path = WF_dir + "/" + shapefiles[wf].name
        WF_lon, WF_lat = Extract_WF_Coordinates(fire_area_path)

        percent = []
        for i in range(len(SD)):
            sd = SD[i]

            # Off info for generators
            out_gen_indices, out_gen_indices_dis = generators_within_safe_distance(gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, sd)

            # Off info for lines
            WF_lines_index, WF_lines_id, WF_line_dis = lines_within_safe_distance(features, WF_lon, WF_lat, sd, weird_coord_structure)
            
            # Find off lines and off transformers
            WF_trans_index = np.intersect1d(WF_lines_index, transformer_index)
            WF_non_trans_index = np.setdiff1d(WF_lines_index, WF_trans_index)

            percent.append((len(WF_non_trans_index), len(WF_trans_index), len(out_gen_indices)))
            print("Finished calculating infomation for percentile: " + str(percentiles[wf]) + " with SD: " + str(SD[i]))
            if display_off_info:
                print("The off information due to the wildfire:")   
                print("Number of off lines: " + str(len(self.WF_non_trans_index)))
                print("Number of off transformers: " + str(len(self.WF_trans_index)))
                print("Number of off generator buses: " + str(len(self.out_gen_indices)) + "\n")    
        
        # Store the data
        data = {
                    str(percentiles[wf]): percent
                }
        
        # Convert the tuples to strings
        data = {key: [str(val) for val in value] for key, value in data.items()}

        # Add new columns for each SD
        new_cols_df = pd.DataFrame(data, index=SD)
        df = pd.concat([df, new_cols_df], axis=1)
    
    if save_csv == True:
        if csv_path == None:
            df.to_csv("./off_info_data.csv", index=True, index_label="Distance / Percentiles")
        else:
            df.to_csv(csv_path, index=True, index_label="Distance / Percentiles")

    return df

# Make plot for the off info from Off_info_vai_SD_Percentile
def Make_off_info_plot_from_CSV(df):
    # Convert cleaned data to tuples
    data = df.map(parse_tuple)

    # Extract dimensions
    rows, cols = data.shape

    # Initialize empty matrices
    A = np.zeros((rows, cols))
    B = np.zeros((rows, cols))
    C = np.zeros((rows, cols))

    # Fill matrices
    for i in range(rows):
        for j in range(cols):
            A[i, j], B[i, j], C[i, j] = data.iloc[i, j]  # Unpacking tuple

    # Compute the fourth matrix
    D = A + B + C

    # Create X, Y meshgrid
    X, Y = np.meshgrid(np.arange(A.shape[1]), np.arange(A.shape[0]))

    # List of matrices and their titles
    matrices = [A, B, C, D]
    titles = ["Off Lines", "Off Transformer", "Off Generators", "Total Off Elements"]

    # # Define custom ticks
    x_ticks = df.columns
    y_ticks = df.index

    # Create a 2×2 subplot layout
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Loop over each rotated matrix and plot it
    for ax, matrix, title in zip(axes.flat, matrices, titles):
        im = ax.imshow(matrix, cmap="viridis", aspect="auto")  # 2D heatmap
        ax.set_title(title)
        ax.set_xlabel("Percentiles")
        ax.set_ylabel("Distances")
        ax.set_xticks(np.linspace(0, matrix.shape[1] - 1, len(x_ticks)))  # Set X ticks
        ax.set_xticklabels(x_ticks)  # Set X tick labels
        ax.set_yticks(np.linspace(0, matrix.shape[0] - 1, len(y_ticks)))  # Set Y ticks
        ax.set_yticklabels(y_ticks)  # Set Y tick labels
        fig.colorbar(im, ax=ax)  # Add colorbar to each plot

    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()


# Collected the information the generators, lines and transformers turned off  to be used in Off_info_plot
def Make_off_info_plot(CATS_lines_dir, CATS_buses_dir, CATS_gens_dir, WF_data_gdfs, WF_lons, WF_lats, out_gen_indixs, WF_lines_indexs, con_names, ylim = None, xlim= None, lon_perturb = .3, lat_perturb = .3):
    
    features, lines_coords, weird_coord_structure = Extract_lines_info(CATS_lines_dir)
    x_segments, y_segments = Lines_segments(lines_coords)

    bus_data = pd.read_csv(CATS_buses_dir) 
    nodes_lat = bus_data.iloc[:, 4].values
    nodes_lon = bus_data.iloc[:, 5].values

    gen_df = pd.read_csv(CATS_gens_dir)  
    all_gen_buses = np.array(gen_df['bus'])
    gen_buses = np.unique(all_gen_buses)

    for k in range(len(WF_data_gdfs)):
        WF_data_gdf = WF_data_gdfs[k]
        WF_lon = WF_lons[k]
        WF_lat = WF_lats[k] 
        out_gen_indix = out_gen_indixs[k]
        WF_lines_index = WF_lines_indexs[k]
        con_name = con_names[k]

        if ylim == None or xlim == None:
            wf_lon_max, wf_lon_min = min(WF_lon), max(WF_lon)
            wf_lat_max, wf_lat_min = min(WF_lat), max(WF_lat)

            xlim = [wf_lon_max - lon_perturb, wf_lon_min + lon_perturb]
            ylim = [wf_lat_max - lat_perturb, wf_lat_min + lat_perturb]

        line_colors = Update_line_colors(lines_coords, WF_lines_index)
        bus_colors = Update_nodes_colors(gen_buses, nodes_lat, nodes_lon, out_gen_indix)

        Off_info_plot(WF_data_gdf, x_segments, y_segments, line_colors, nodes_lon, nodes_lat, bus_colors, ylim = ylim, xlim= xlim, title = con_name, loc='best')

# Create a set containing vectors with the information of the lines coordinates for each line.
def Lines_segments(lines_coords):
    # Initialize lists for segments
    x_segments = []
    y_segments = []

    for i, line_coords in enumerate(lines_coords):
        num_coords = len(line_coords)
        x_segment = []
        y_segment = []
        
        for j in range(num_coords - 1):
            lon_src, lat_src = line_coords[j]
            lon_dest, lat_dest = line_coords[j + 1]
            x_segment.extend([lon_src, lon_dest])
            y_segment.extend([lat_src, lat_dest])
        
        x_segments.append(x_segment)
        y_segments.append(y_segment)

    return x_segments, y_segments

# Update/Define the colors for each line
def Update_line_colors(lines_coords, WF_lines_index):
    line_colors = []
    for i, line_coords in enumerate(lines_coords):
        num_coords = len(line_coords)
        
        # Check what lines are within safe distance of the WF
        if i in WF_lines_index:
            line_colors.append("red")
        else:
            line_colors.append("green")

    return line_colors

# Update/Define the colors for each node/bus
def Update_nodes_colors(gen_buses, nodes_lat, nodes_lon, out_gen_indices):
    bus_colors = []
    load_bus_color = "black"
    gen_bus_color = "cyan"
    gen_bus_WF_color = "red"

    if len(out_gen_indices) > 0:
        for i in range(len(nodes_lat)):
            j = i + 1
            if j in gen_buses[out_gen_indices]:
                bus_colors.append(gen_bus_WF_color)
            elif j in gen_buses:
                bus_colors.append(gen_bus_color)
            else:
                bus_colors.append(load_bus_color)
    else:
        for i in range(len(nodes_lat)):
            j = i + 1
            if j in gen_buses:
                bus_colors.append(gen_bus_color)
            else:
                bus_colors.append(load_bus_color)

    return np.array(bus_colors)

# Make plot showing the generators, lines and transformers turned off on the map of California
def Off_info_plot(fire_area_gdf, x_segments, y_segments, line_colors, nodes_lon, nodes_lat, bus_colors, ylim = None, xlim= None, title = None, loc='best', node_size=10):
    lines_legends = ['On lines', 'Off lines (WF)']
    lines_colors_legends = [ "green", "red"]

    nodes_legends = ['loads', 'On generators', 'Off generators (WF)']
    nodes_colors_legends = ["black", "cyan", "red"]

    # Plot the map
    fig, ax = plt.subplots(figsize=(10, 10))
    # Set axis bounds

    if xlim != None and ylim != None: 
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

    fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326)  
    fire_area_4326_gdf.plot(ax=ax, color="gray", zorder=1)

    for x, y, color in zip(x_segments, y_segments, line_colors):
        ax.plot(x, y, color=color, linewidth=1, zorder=2)

    # Overlay the nodes with a higher zorder
    ax.scatter(nodes_lon, nodes_lat, color=bus_colors, s=node_size, zorder=3)

    # Create legend handles for nodes
    node_legend_handles = [mlines.Line2D([], [], color=color, marker='o', linestyle='None', markersize=8, label=name) 
                        for name, color in zip(nodes_legends, nodes_colors_legends)]

    # Create legend handles for edges
    edge_legend_handles = [mlines.Line2D([], [], color=color, linewidth=2, label=name) 
                        for name, color in zip(lines_legends, lines_colors_legends)]

    # Add legend
    ax.legend(handles=node_legend_handles + edge_legend_handles, loc=loc)

    ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string(), zorder=0)

    if title != None:
        ax.set_title(title)

    # Customize the plot
    plt.show()


# Functions to create multiple figures depending on varying SD and percentile or just SD
# Made for probabilistic wildfire data
def Make_figs_elm_near_WF_vary_Percentiles_SD(percentiles, SD, WF_dir, shapefiles, x_segments, y_segments, line_colors, nodes_lon, nodes_lat, bus_colors, ylim = None, xlim= None, loc='best', node_size=10, save_fig = False, fig_dir = None, lon_perturb = .3, lat_perturb = .3):
    lines_legends = ['On lines', 'Off lines (WF)']
    lines_colors_legends = [ "green", "red"]

    nodes_legends = ['loads', 'On generators', 'Off generators (WF)']
    nodes_colors_legends = ["black", "cyan", "red"]
    for i in range(len(percentiles)):
        fire_area_path = WF_dir + "/" + shapefiles[i].name
        fire_area_gdf = gpd.read_file(fire_area_path)
        fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326) 

        if ylim == None or xlim == None:
            WF_lon, WF_lat = Extract_WF_Coordinates(fire_area_path)
            wf_lon_max, wf_lon_min = min(WF_lon), max(WF_lon)
            wf_lat_max, wf_lat_min = min(WF_lat), max(WF_lat)

            xlim = [wf_lon_max - lon_perturb, wf_lon_min + lon_perturb]
            ylim = [wf_lat_max - lat_perturb, wf_lat_min + lat_perturb]
            
        for sd in range(len(SD)):
            # Create the plot
            fig, ax = plt.subplots(figsize=(10, 10))

            # Set axis bounds
            ax.set_ylim(ylim)
            ax.set_xlim(xlim)
 
            fire_area_4326_gdf.plot(ax=ax, color="gray", zorder=1)

            # Plot the segments
            for x, y, color in zip(x_segments, y_segments, line_colors[sd + (i - 1) * len(SD)]):
                ax.plot(x, y, color=color, linewidth=1, zorder=2)

            # Overlay the nodes with a higher zorder
            ax.scatter(nodes_lon, nodes_lat, color=bus_colors[sd + (i - 1) * len(SD)], s=10, zorder=3)

            # Create legend handles for nodes
            node_legend_handles = [mlines.Line2D([], [], color=color, marker='o', linestyle='None', markersize=8, label=name) 
                                for name, color in zip(nodes_legends, nodes_colors_legends)]

            # Create legend handles for edges
            edge_legend_handles = [mlines.Line2D([], [], color=color, linewidth=2, label=name) 
                                for name, color in zip(lines_legends, lines_colors_legends)]

            # Add legend
            ax.legend(handles=node_legend_handles + edge_legend_handles, loc='upper right')

            ax.set_title("Percentile: " + str(percentiles[i]) + ", SD: " + str(SD[sd]))

            ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string(), zorder=0)

            if save_fig:
                if fig_dir == None:
                    os.makedirs("plot_dir", exist_ok=True)
                    plt.savefig("./plot_dir/Percentile: " + str(percentiles[i]) + ", SD: " + str(SD[sd]), dpi=300, bbox_inches='tight') 
                else:
                    os.makedirs(fig_dir, exist_ok=True)
                    plt.savefig("./" + fig_dir + "/Percentile: " + str(percentiles[i]) + ", SD: " + str(SD[sd]), dpi=300, bbox_inches='tight') 

             # Display the plot
            plt.show()

def Make_Series_of_figs_near_WF_vary_Percentiles_SD(percentiles, SD, WF_dir, CATS_buses_dir, CATS_gens_dir, CATS_lines_dir, ylim = None, xlim= None, loc='best', node_size=10, save_fig = False, fig_dir = None, lon_perturb = .3, lat_perturb = .3):   
    if not isinstance(percentiles, list): 
        percentiles = [percentiles]

    if not isinstance(SD, list): 
        SD = [SD]
    
    bus_df = pd.read_csv(CATS_buses_dir)  

    nodes_lat = bus_df.iloc[:, 4].values
    nodes_lon = bus_df.iloc[:, 5].values

    gen_buses, gen_buses_units, gen_buses_lon, gen_buses_lat = Extract_gen_info(CATS_gens_dir)

    features, lines_coords, weird_coord_structure = Extract_lines_info(CATS_lines_dir)

    shapefiles = ID_WF_files(WF_dir, percentiles)

    if len(shapefiles) != len(percentiles):
        print("Warning: One of the percentiles data was not within the directory")

    bus_colors = []
    line_colors = []

    x_segments, y_segments = Lines_segments(lines_coords)

    for wf in range(len(shapefiles)):
        fire_area_path = WF_dir + "/" + shapefiles[wf].name
        WF_lon, WF_lat = Extract_WF_Coordinates(fire_area_path)

        for i in range(len(SD)):
            sd = SD[i]

            out_gen_indices, out_gen_indices_dis = generators_within_safe_distance(gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, sd)
            
            WF_lines_index, WF_lines_id, WF_line_dis = lines_within_safe_distance(features, WF_lon, WF_lat, sd, weird_coord_structure)

            bus_colors.append(Update_nodes_colors(gen_buses, nodes_lat, nodes_lon, out_gen_indices))

            color = Update_line_colors(lines_coords, WF_lines_index)

            line_colors.append(color)

    Make_figs_elm_near_WF_vary_Percentiles_SD(percentiles, SD, WF_dir, shapefiles, x_segments, y_segments, line_colors, nodes_lon, nodes_lat, bus_colors, ylim = ylim, xlim= xlim, loc=loc, node_size=node_size, save_fig = save_fig, fig_dir = fig_dir, lon_perturb = lon_perturb, lat_perturb = lat_perturb)

# Made for historic wildfire data
def Make_figs_elm_near_WF_vary_SD(SD, WF_dir, x_segments, y_segments, line_colors, nodes_lon, nodes_lat, bus_colors, ylim = None, xlim= None, loc='best', node_size=10, save_fig = False, fig_dir = None, lon_perturb = .3, lat_perturb = .3):
    lines_legends = ['On lines', 'Off lines (WF)']
    lines_colors_legends = [ "green", "red"]

    nodes_legends = ['loads', 'On generators', 'Off generators (WF)']
    nodes_colors_legends = ["black", "cyan", "red"]
        # Create the plot

    fire_area_gdf = gpd.read_file(WF_dir)
    fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326)  

    WF_lon, WF_lat = Extract_WF_Coordinates(WF_dir)
    if ylim == None or xlim == None:
        wf_lon_max, wf_lon_min = min(WF_lon), max(WF_lon)
        wf_lat_max, wf_lat_min = min(WF_lat), max(WF_lat)

        xlim = [wf_lon_max - lon_perturb, wf_lon_min + lon_perturb]
        ylim = [wf_lat_max - lat_perturb, wf_lat_min + lat_perturb]

    for sd in range(len(SD)):
        fig, ax = plt.subplots(figsize=(10, 10))

        fire_area_4326_gdf.plot(ax=ax, color="gray", zorder=1)
        
        # Set axis bounds
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)

        # Plot the segments
        for x, y, color in zip(x_segments, y_segments, line_colors[sd ]):
            ax.plot(x, y, color=color, linewidth=1, zorder=2)

        # Overlay the nodes with a higher zorder
        ax.scatter(nodes_lon, nodes_lat, color=bus_colors[sd ], s=10, zorder=3)

        # Create legend handles for nodes
        node_legend_handles = [mlines.Line2D([], [], color=color, marker='o', linestyle='None', markersize=8, label=name) 
                            for name, color in zip(nodes_legends, nodes_colors_legends)]

        # Create legend handles for edges
        edge_legend_handles = [mlines.Line2D([], [], color=color, linewidth=2, label=name) 
                            for name, color in zip(lines_legends, lines_colors_legends)]

        # Add legend
        ax.legend(handles=node_legend_handles + edge_legend_handles, loc='upper right')

        ax.set_title("SD: " + str(SD[sd]))

        ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string(), zorder=0)

        if save_fig:
            if fig_dir == None:
                os.makedirs("plot_dir", exist_ok=True)
                plt.savefig("./plot_dir/SD_" + str(SD[sd]), dpi=300, bbox_inches='tight') 
            else:
                os.makedirs(fig_dir, exist_ok=True)
                plt.savefig("./" + fig_dir + "SD_" + str(SD[sd]), dpi=300, bbox_inches='tight')             

        # Display the plot
        plt.show()

def Make_Series_of_figs_near_WF_vary_SD(SD, WF_dir, CATS_buses_dir, CATS_gens_dir, CATS_lines_dir, ylim = None, xlim= None, loc='best', node_size=10, save_fig = False, fig_dir = None, lon_perturb = .3, lat_perturb = .3):   
    
    bus_df = pd.read_csv(CATS_buses_dir)  

    nodes_lat = bus_df.iloc[:, 4].values
    nodes_lon = bus_df.iloc[:, 5].values

    gen_buses, gen_buses_units, gen_buses_lon, gen_buses_lat = Extract_gen_info(CATS_gens_dir)

    features, lines_coords, weird_coord_structure = Extract_lines_info(CATS_lines_dir)

    fire_area_gdf = gpd.read_file(WF_dir)
    WF_lon, WF_lat = Extract_WF_Coordinates(WF_dir)

    bus_colors = []
    line_colors = []

    x_segments, y_segments = Lines_segments(lines_coords)

    for i in range(len(SD)):
        sd = SD[i]

        out_gen_indices, out_gen_indices_dis = generators_within_safe_distance(gen_buses_lon, gen_buses_lat, WF_lon, WF_lat, sd)
        
        WF_lines_index, WF_lines_id, WF_line_dis = lines_within_safe_distance(features, WF_lon, WF_lat, sd, weird_coord_structure)

        bus_colors.append(Update_nodes_colors(gen_buses, nodes_lat, nodes_lon, out_gen_indices))

        color = Update_line_colors(lines_coords, WF_lines_index)

        line_colors.append(color)

    Make_figs_elm_near_WF_vary_SD(SD, WF_dir, x_segments, y_segments, line_colors, nodes_lon, nodes_lat, bus_colors, ylim = ylim, xlim= xlim, loc=loc, node_size=node_size, save_fig = save_fig, fig_dir = fig_dir, lon_perturb = lon_perturb, lat_perturb = lat_perturb)


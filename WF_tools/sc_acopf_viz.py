from required_funtions import *
from Make_fig_funtions import *
from parse_data import *
from validation_functions import *

# Classes use to visualize the different data provide by ExaJuGo
class WF_figures:
    def __init__(self, CATS_lines_circuits_dir, CATS_transformer_dir, CATS_gens_dir, CATS_lines_dir, CATS_buses_dir,
                    CATS_case_dir):
        """Initialize the class with variables."""
        self.CATS_lines_circuits_dir = CATS_lines_circuits_dir
        self.CATS_transformer_dir = CATS_transformer_dir
        self.CATS_gens_dir = CATS_gens_dir
        self.CATS_lines_dir = CATS_lines_dir
        self.CATS_buses_dir = CATS_buses_dir
        self.CATS_case_dir = CATS_case_dir

        self.transformer_index                                                       = ID_transformer(CATS_transformer_dir)
        self.features, self.lines_coords, self.weird_coord_structure                 = Extract_lines_info(CATS_lines_dir)
        self.bus_info, self.nodes_lat, self.nodes_lon                                = ID_bus_info(CATS_buses_dir, self.features)
        self.gen_buses, self.gen_buses_units, self.gen_buses_lon, self.gen_buses_lat = Extract_gen_info(CATS_gens_dir)
        self.from_bus, self.to_bus                                                   = endpoint_info(self.features)
        self.bus_data                                                                = pd.read_csv(CATS_buses_dir) 
        self.nodes_lat                                                               = self.bus_data.iloc[:, 4].values
        self.nodes_lon                                                               = self.bus_data.iloc[:, 5].values
        self.Define_bus_info()

        with open(CATS_case_dir, "r") as f:
            self.cats_case = json.load(f)

        self.WF_lon = []
        self.WF_lat = []
        self.WF_lines_index = []

        self.x_segments = []
        self.y_segments = []
        self.line_colors = []

        self.EF_base_pf = []
        self.EF_contingencies_pf = []
        self.ACOPF_base = []
        self.avg_EF_pf_bc = []
        self.avg_EF_pf_con = []
        self.percent_delta_pf = []
        self.delta = []
        
        self.EF_slacks_data = []
        self.line_slacks = []

        self.off_nodes_index = []
        self.off_nodes = []
        self.off_gen_nodes = [] 
        self.lost_lines_index = []
        self.num_lost_lines = []

        self.A_Power_demand = []
        self.R_Power_demand = []

        self.line_safe_color = "none"
        self.line_danger_color = "grey"

    # Create bus object containt bus id, index of the lines connected to the bus, total number of lines connected to the bus
    def Define_bus_info(self):
        bus_id = np.array(self.bus_data['bus_i'])

        bus_lines_index = []
        total_bus_lines_index = []
        for x in bus_id:
            # Get indices where x appears in both arrays
            indices1 = np.where(self.from_bus == x)[0]  
            indices2 = np.where(self.to_bus == x)[0]  

            # Combine indices into a single 1D array
            combined_indices = np.concatenate((indices1, indices2))
            bus_lines_index.append(combined_indices)
            total_bus_lines_index.append(len(combined_indices))

        self.bus_info = {'id':bus_id, 'bus_lines_index': bus_lines_index, 'total_lines': np.array(total_bus_lines_index)}

    # Define the longitude and latitude cooridinates of the wildfire of interest
    def Update_WF_coords(self, WF_lon, WF_lat):
        self.WF_lon = WF_lon
        self.WF_lat = WF_lat

    # Find the lines, transformer and generators turned off in the SD mile radius
    def Create_WF_lines_index(self, sd):
        if len(self.WF_lon) == 0:
            raise ValueError("WF_lon and WF_lat was not prodived. Use the function Update_WF_coords to provide them")

        self.WF_lines_index, _, _ = lines_within_safe_distance(self.features, self.WF_lon, self.WF_lat, sd, self.weird_coord_structure)
        self.out_gen_indices, _ = generators_within_safe_distance(self.gen_buses_lon, self.gen_buses_lat, self.WF_lon, self.WF_lat, sd)
        self.out_node_indices, _ = generators_within_safe_distance(self.nodes_lon, self.nodes_lat, self.WF_lon, self.WF_lat, sd)

        # ID the indices that belong to transformer or lines (non_trans)
        self.WF_trans_index = np.intersect1d(self.WF_lines_index, self.transformer_index)
        self.WF_non_trans_index = np.setdiff1d(self.WF_lines_index, self.WF_trans_index)

        # print number of off info
        print("The off information due to the wildfire:")
        print("Number of off lines: " + str(len(self.WF_non_trans_index)))
        print("Number of off transformers: " + str(len(self.WF_trans_index)))
        print("Number of off generator buses: " + str(len(self.out_gen_indices)) + "\n")

    # Updated the colors associated to the off line due to the wildfire and the on lines
    def Update_line_colors(self, safe_color, danger_color):
        self.line_colors = []
        for i, line_coords in enumerate(self.lines_coords):
            num_coords = len(line_coords)
            
            # Check what lines are within safe distance of the WF
            if i in self.WF_lines_index:
                self.line_colors.append(danger_color)
            else:
                self.line_colors.append(safe_color)

    # Create a set containing vectors with the information of the lines coordinates for each line.
    def Lines_segments(self):
        # Initialize lists for segments
        self.x_segments = []
        self.y_segments = []
        self.line_colors = []

        safe_color = self.line_safe_color
        danger_color = self.line_danger_color

        for i, line_coords in enumerate(self.lines_coords):
            num_coords = len(line_coords)
            x_segment = []
            y_segment = []
            
            # Check what lines are within safe distance of the WF
            if i in self.WF_lines_index:
                for j in range(num_coords - 1):
                    lon_src, lat_src = line_coords[j]
                    lon_dest, lat_dest = line_coords[j + 1]
                    x_segment.extend([lon_src, lon_dest])
                    y_segment.extend([lat_src, lat_dest])
                self.line_colors.append(danger_color)
            else:
                for j in range(num_coords - 1):
                    lon_src, lat_src = line_coords[j]
                    lon_dest, lat_dest = line_coords[j + 1]
                    x_segment.extend([lon_src, lon_dest])
                    y_segment.extend([lat_src, lat_dest])
                self.line_colors.append(safe_color)
            
            self.x_segments.append(x_segment)
            self.y_segments.append(y_segment)

    # Save the power flow data in the object
    def Create_Power_Flow_data(self, EF_filename, small_value = 1e-5):
        # EF_base_pf: The power flow data for the basecase
        # EF_contingencies_pf: The power flow data after a contingency
        self.EF_base_pf, self.EF_contingencies_pf = parse_power_data(EF_filename)

        # The average of the power going oppisite direction in the same lines
        self.avg_EF_pf_bc  = (np.abs(self.EF_base_pf['p_from']) + np.abs(self.EF_base_pf['p_to']) ) *.5
        self.avg_EF_pf_con = (np.abs(self.EF_contingencies_pf[0]['p_from']) + np.abs(self.EF_contingencies_pf[0]['p_to']) ) *.5

        # The difference between the power flows
        delta = self.avg_EF_pf_con - self.avg_EF_pf_bc
        self.percent_delta_pf = 0 * self.avg_EF_pf_con

        # Percentage of the power flowcompared the max between the basecase and contingency power flow
        for i in range(len(delta)):
            self.percent_delta_pf[i] = 100 * delta[i] / max(abs(self.avg_EF_pf_con[i]), abs(self.avg_EF_pf_bc[i]), small_value)

        self.delta = delta

    # Plot power flow data with color bar
    def Make_Delta_PF_plot(self, fire_area_gdf, ylim = None, xlim= None, shrink = 1.0, title = None, loc='best'):
        fig, ax = plt.subplots(figsize=(8, 6))

        # Set axis bounds
        if xlim != None and ylim != None:
            ax.set_ylim(ylim)
            ax.set_xlim(xlim)

        # legends for the edges
        edge_legend_handles = mlines.Line2D([], [], color="dimgrey", linewidth=2, label="lines turned off") 

        # the wildfire panda date frame
        fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326)
        fire_area_4326_gdf.plot(ax=ax, color="gray", zorder=1)

        # Define colormap and normalization
        vmin, vmax  = min(self.delta), max(self.delta)

        # Compute the normalized position of zero
        zero_position = (0 - vmin) / (vmax - vmin)

        # create custom color bar
        cmap = mcolors.LinearSegmentedColormap.from_list(
            "custom_cmap", 
            [(0.0, "blue"), (zero_position, "yellow"), (1.0, "red")]
        )

        # Normaliza the color bar so data lie between the min and max
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # Plot each path as a polyline
        perm = np.argsort(abs(self.delta))
        for i in perm:
            # Each x and y here are full paths, not just start and end points
            x_path = self.x_segments[i]
            y_path = self.y_segments[i]
            if i in self.WF_lines_index:
                color = "dimgrey"
            else:
                value = self.delta[i]  # The corresponding value for color mapping
                normalized_value = norm(value)  # Convert to [0,1] range
                color = cmap(normalized_value)  # Get the correct color
            
            # Plot the path (polyline) with the corresponding color
            ax.plot(x_path, y_path, color=color, linewidth=0.5, zorder=2)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])  # Required for colorbar
        cbar = plt.colorbar(sm, ax=ax, shrink = shrink, aspect=30)
        cbar.set_label("Delta Active Power Flow (MW)")

        if title != None:
            ax.set_title(title)

        ax.legend(handles= [edge_legend_handles], loc=loc)

        ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string(), zorder=0)

        plt.show()

    # Same as Make_Delta_PF_plot but plot the power flow data as percent
    def Make_Delta_PF_percent_plot(self, fire_area_gdf, ylim = None, xlim= None, shrink = 1.0, title = None, loc='best'):
        fig, ax = plt.subplots(figsize=(8, 6))

        # Set axis bounds
        if xlim != None and ylim != None:
            ax.set_ylim(ylim)
            ax.set_xlim(xlim)

        edge_legend_handles = mlines.Line2D([], [], color="dimgrey", linewidth=2, label="lines turned off") 

        fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326)
        fire_area_4326_gdf.plot(ax=ax, color="gray", zorder=1)

        # Define colormap and normalization
        vmin, vmax = -100, 100

        # Compute the normalized position of zero
        zero_position = (0 - vmin) / (vmax - vmin)

        cmap = mcolors.LinearSegmentedColormap.from_list(
            "custom_cmap", 
            [(0.0, "blue"), (zero_position, "yellow"), (1.0, "red")]
        )

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # Plot each path as a polyline
        perm = np.argsort(abs(self.percent_delta_pf))
        for i in perm:
            # Each x and y here are full paths, not just start and end points
            x_path = self.x_segments[i]
            y_path = self.y_segments[i]
            if i in self.WF_lines_index:
                color = "dimgrey"
            else:
                value = self.percent_delta_pf[i]  # The corresponding value for color mapping
                normalized_value = norm(value)  # Convert to [0,1] range
                color = cmap(normalized_value)  # Get the correct color
            
            # Plot the path (polyline) with the corresponding color
            ax.plot(x_path, y_path, color=color, linewidth=0.5, zorder=2)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])  # Required for colorbar
        cbar = plt.colorbar(sm, ax=ax, shrink = shrink, aspect=30)
        cbar.set_label("Delta Active Power Flow in Percent")

        if title != None:
            ax.set_title(title)

        ax.legend(handles= [edge_legend_handles], loc=loc)

        ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string(), zorder=0)

        plt.show()

    # Save the slack data to the object
    def Define_slack_data(self, SCACOPF_slack_filename):
        self.EF_slacks_data = collect_EF_slack_data(SCACOPF_slack_filename)
        self.line_slacks = np.concatenate([np.array(self.EF_slacks_data['contingency'][0]['line']['slack']) , np.array(self.EF_slacks_data['contingency'][0]['transformer']['slack'])])
        self.slack_buses_index = np.where(self.EF_slacks_data['contingency'][0]['node']['slack'] > 0 )[0]

    # Plot power flow slacks from ExaJuGo
    def Make_Slack_PF_plot(self, fire_area_gdf, ylim = None, xlim= None, shrink = 1.0, title = None, loc='best'):
        fig, ax = plt.subplots(figsize=(8, 6))

        # Set axis bounds
        if xlim != None and ylim != None:
            ax.set_ylim(ylim)
            ax.set_xlim(xlim)

        # legends for the edges
        edge_legend_handles = mlines.Line2D([], [], color="grey", linewidth=2, label="off lines") 

        # the wildfire panda date frame
        fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326)
        fire_area_4326_gdf.plot(ax=ax, color="gray", zorder=1)

        # take the min and max value if they are nonzero and set them to 0 and 1
        if max(self.line_slacks) == 0:
            vmin, vmax = 0, 1
        else:
            vmin, vmax = 0, max(self.line_slacks)

        # create custom color bar
        cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["yellow", "red"])

        # Normaliza the color bar so data lie between the min and max
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # Plot each path as a polyline
        for i in range(len(self.x_segments)):
            # Each x and y here are full paths, not just start and end points
            x_path = self.x_segments[i]
            y_path = self.y_segments[i]
            if i in self.WF_lines_index:
                color = "grey"
            else:
                value = self.line_slacks[i]  # The corresponding value for color mapping
                normalized_value = norm(value)  # Convert to [0,1] range
                color = cmap(normalized_value)  # Get the correct color
            
            # Plot the path (polyline) with the corresponding color
            ax.plot(x_path, y_path, color=color, linewidth=0.5, zorder=2)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])  # Required for colorbar
        cbar = plt.colorbar(sm, ax=ax, shrink = shrink, aspect=30)
        cbar.set_label("Slack")

        if title != None:
            ax.set_title(title)

        ax.legend(handles= [edge_legend_handles], loc=loc)

        ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string(), zorder=0)

        plt.show()

    # Save the information of off elements
    def ID_off_info(self):
        self.off_nodes_index, self.off_nodes, self.off_gen_nodes, self.lost_lines_index, self.num_lost_lines = ID_off_and_lost_info(self.bus_info, self.WF_lines_index, self.gen_buses)

    # ID the different type of buses
    # gen_buses_only_index: Index of Buses containing only generators
    # gen_load_buses_index: Index of Buses containing generators and loads
    # neither_gen_load_buses_index: Index of Buses contaning no generators or loads
    # load_buses_only_index: Index of Buses containing only loads
    # Slack_: ID which of those buses have slack
    def Node_type_Index(self):
        A_Power_demand = []
        R_Power_demand = []

        # The power demand from each bus
        for bus in self.cats_case['bus']:
            A_Power_demand.append(bus['Pd'])
            R_Power_demand.append(bus['Qd'])

        # The load buses and their index
        self.load_buses_index = np.union1d(np.where(np.array(A_Power_demand) > 0)[0], np.where(np.array(R_Power_demand) > 0)[0])
        self.load_buses = self.load_buses_index + 1

        self.gen_load_buses = np.intersect1d(self.gen_buses, self.load_buses)

        self.neither_gen_load_buses = np.setdiff1d( self.bus_data['bus_i'],np.union1d(self.gen_buses, self.load_buses))
        self.load_buses_only = np.setdiff1d(self.load_buses, self.gen_buses)
        self.gen_buses_only = np.setdiff1d(self.gen_buses, self.load_buses)

        # Assuming the buses enumeration start at 1 so index is shifting the enumeration
        self.gen_buses_only_index = self.gen_buses_only - 1
        self.gen_load_buses_index = self.gen_load_buses - 1
        self.neither_gen_load_buses_index = self.neither_gen_load_buses - 1
        self.load_buses_only_index = self.load_buses_only - 1

        # ID which buses have slack
        self.slack_gen_load_buses_index = np.intersect1d(self.EF_slacks_data['basecase']['node']['i'][self.slack_buses_index] - 1, self.gen_load_buses_index)
        self.slack_neither_gen_load_buses_index = np.intersect1d(self.EF_slacks_data['basecase']['node']['i'][self.slack_buses_index] - 1, self.neither_gen_load_buses_index)
        self.slack_gen_buses_only_index = np.intersect1d(self.EF_slacks_data['basecase']['node']['i'][self.slack_buses_index] - 1, self.gen_buses_only_index)
        self.slack_load_buses_only_index = np.intersect1d(self.EF_slacks_data['basecase']['node']['i'][self.slack_buses_index] - 1, self.load_buses_only_index)

    # Update the color of the node for plot
    def Update_nodes_colors(self, load_bus_color = "black", gen_bus_color = "cyan", gen_bus_WF_color = "red"):
        bus_colors = []

        if len(self.out_gen_indices) > 0:
            for i in range(len(self.nodes_lat)):
                j = i + 1
                if j in self.gen_buses[self.out_gen_indices]:
                    bus_colors.append(gen_bus_WF_color)
                elif j in self.gen_buses:
                    bus_colors.append(gen_bus_color)
                else:
                    bus_colors.append(load_bus_color)
        else:
            for i in range(len(self.nodes_lat)):
                j = i + 1
                if j in self.gen_buses:
                    bus_colors.append(gen_bus_color)
                else:
                    bus_colors.append(load_bus_color)

        self.bus_colors = np.array(bus_colors)

    # Find nodes affected by the wildfire and given their associated color
    def Node_affected(self):
        bus_colors = []
        affected = "red"
        not_affected = "green"

        if len(self.out_gen_indices) != 0:
            affected_buses = np.union1d(self.gen_buses[self.out_gen_indices], np.union1d(self.off_nodes, self.lost_lines_index + 1))
        else:
            affected_buses = np.union1d(self.off_nodes, self.lost_lines_index + 1)

        for i in range(len(self.bus_info['id'])):
            j = i + 1
            if j in affected_buses:
                bus_colors.append(affected)
            else:
                bus_colors.append(not_affected)

        self.bus_colors = np.array(bus_colors)

    # Create the legends for the lines
    def Line_lengends(self, safe_color, danger_color):
        self.Update_line_colors(safe_color, danger_color)

        # Check if line is on or it was turned off due to wildfire
        possible_lines_legends = ['On lines', 'Off lines (WF)']
        possible_lines_colors_legends = [safe_color, danger_color]

        self.lines_legends = []
        self.lines_colors_legends = []

        for lcl in range(len(possible_lines_colors_legends)):
            if len(np.where(np.array(self.line_colors) == possible_lines_colors_legends[lcl])[0]) > 0 and possible_lines_colors_legends[lcl] != 'none':
                self.lines_legends.append(possible_lines_legends[lcl])
                self.lines_colors_legends.append(possible_lines_colors_legends[lcl])

    # Make the plots of the nodes with slacks
    def Make_Node_Slack_Plot(self, only_colorbar, fire_area_gdf, ylim = None, xlim= None, 
                                shrink = 1.0, title = None, loc='best'):
        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 10))

        # Color map with colors from Wistia
        cmap = plt.cm.Wistia  

        # Wildfire geopanda object
        fire_area_4326_gdf = fire_area_gdf.to_crs(epsg=4326)  
        fire_area_4326_gdf.plot(ax=ax, color="gray", zorder=1)

        # Plot the segments
        for x, y, color in zip(self.x_segments, self.y_segments, self.line_colors):
            ax.plot(x, y, color=color, linewidth=.5, zorder=2)

        # Create legend handles for edges
        edge_legend_handles = [mlines.Line2D([], [], color=color, linewidth=2, label=name) 
                            for name, color in zip(self.lines_legends, self.lines_colors_legends)]

        # If we only want to see the buses with slacks or all of the nodes
        if only_colorbar:
            if xlim != None and ylim != None:
                ax.set_ylim(ylim)
                ax.set_xlim(xlim)

            # Scatter numerical nodes with the custom colormap
            num_no_load_gen_scatter = ax.scatter(
                [self.nodes_lon[n] for n in self.slack_neither_gen_load_buses_index],
                [self.nodes_lat[n] for n in self.slack_neither_gen_load_buses_index],
                c=[n for n in self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_neither_gen_load_buses_index]],
                cmap=cmap, s=40, marker='X', linewidths=1, zorder=3
            )
            num_load_gen_scatter = ax.scatter(
                [self.nodes_lon[n] for n in self.slack_gen_load_buses_index],
                [self.nodes_lat[n] for n in self.slack_gen_load_buses_index],
                c=[n for n in self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_gen_load_buses_index]],
                cmap=cmap, s=40, marker='*', zorder=3
            )
            num_load_scatter = ax.scatter(
                [self.nodes_lon[n] for n in self.slack_load_buses_only_index],
                [self.nodes_lat[n] for n in self.slack_load_buses_only_index],
                c=[n for n in self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_load_buses_only_index]],
                cmap=cmap, s=15, marker='^', linewidths=2, zorder=3
            )
            num_gen_scatter = ax.scatter(
                [self.nodes_lon[n] for n in self.slack_gen_buses_only_index],
                [self.nodes_lat[n] for n in self.slack_gen_buses_only_index],
                c=[n for n in self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_gen_buses_only_index]],
                cmap=cmap, s=20, marker='o', zorder=3
            )
            no_load_no_gen_marker_legend = mlines.Line2D([], [], color='black', marker='X', linestyle='None', markersize=8, label="no load and no gen buses", markerfacecolor='none', markeredgecolor='black')
            load_gen_marker_legend = mlines.Line2D([], [], color='black', marker='*', linestyle='None', markersize=8, label="load and gen buses", markerfacecolor='none', markeredgecolor='black')
            load_marker_legend = mlines.Line2D([], [], color='black', marker='^', linestyle='None', markersize=8, label="load buses", markerfacecolor='none', markeredgecolor='black')
            gen_marker_legend = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=8, label="gen buses", markerfacecolor='none', markeredgecolor='black')
            # Add legend
            ax.legend(handles= [gen_marker_legend] + [load_marker_legend] + [load_gen_marker_legend] + [no_load_no_gen_marker_legend] + edge_legend_handles, loc=loc)

        else:
            if xlim != None and ylim != None:
                ax.set_ylim(ylim)
                ax.set_xlim(xlim)

            ax.scatter(
                [self.nodes_lon[n] for n in self.neither_gen_load_buses_index],
                [self.nodes_lat[n] for n in self.neither_gen_load_buses_index],
                color=[self.bus_colors[n] for n in self.neither_gen_load_buses_index], s=40, marker='X', linewidths=1, zorder=5
            )
            ax.scatter(
                [self.nodes_lon[n] for n in self.gen_load_buses_index],
                [self.nodes_lat[n] for n in self.gen_load_buses_index],
                color=[self.bus_colors[n] for n in self.gen_load_buses_index], s=40, marker='*', linewidths=2, zorder=5
            )
            ax.scatter(
                [self.nodes_lon[n] for n in self.load_buses_only_index],
                [self.nodes_lat[n] for n in self.load_buses_only_index],
                color=[self.bus_colors[n] for n in self.load_buses_only_index], s=15, marker='^', linewidths=2, zorder=3
            )
            ax.scatter(
                [self.nodes_lon[n] for n in self.gen_buses_only_index],
                [self.nodes_lat[n] for n in self.gen_buses_only_index],
                color=[self.bus_colors[n] for n in self.gen_buses_only_index], s=20, marker='o', linewidths=2, zorder=4
            )

            # Scatter numerical nodes with the custom colormap
            num_no_load_gen_scatter = ax.scatter(
                [self.nodes_lon[n] for n in self.slack_neither_gen_load_buses_index],
                [self.nodes_lat[n] for n in self.slack_neither_gen_load_buses_index],
                c=[n for n in self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_neither_gen_load_buses_index]],
                cmap=cmap, s=40, marker='X', linewidths=1, zorder=6
            )
            num_load_gen_scatter = ax.scatter(
                [self.nodes_lon[n] for n in self.slack_gen_load_buses_index],
                [self.nodes_lat[n] for n in self.slack_gen_load_buses_index],
                c=[n for n in self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_gen_load_buses_index]],
                cmap=cmap, s=40, marker='*', zorder=6
            )
            num_load_scatter = ax.scatter(
                [self.nodes_lon[n] for n in self.slack_load_buses_only_index],
                [self.nodes_lat[n] for n in self.slack_load_buses_only_index],
                c=[n for n in self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_load_buses_only_index]],
                cmap=cmap, s=15, marker='^', linewidths=2, zorder=6
            )
            num_gen_scatter = ax.scatter(
                [self.nodes_lon[n] for n in self.slack_gen_buses_only_index],
                [self.nodes_lat[n] for n in self.slack_gen_buses_only_index],
                c=[n for n in self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_gen_buses_only_index]],
                cmap=cmap, s=20, marker='o', zorder=6
            )
            no_load_no_gen_marker_legend = mlines.Line2D([], [], color='black', marker='X', linestyle='None', markersize=8, label="no load or gen buses", markerfacecolor='none', markeredgecolor='black')
            load_gen_marker_legend = mlines.Line2D([], [], color='black', marker='*', linestyle='None', markersize=8, label="load and gen buses", markerfacecolor='none', markeredgecolor='black')
            load_marker_legend = mlines.Line2D([], [], color='black', marker='^', linestyle='None', markersize=8, label="load buses", markerfacecolor='none', markeredgecolor='black')
            gen_marker_legend = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=8, label="gen buses", markerfacecolor='none', markeredgecolor='black')
            # Create legend handles for nodes
            affected_legend = mlines.Line2D([], [], color='red', marker='o', linestyle='None', markersize=8, label="Affected by WF") 
            not_affected_legend = mlines.Line2D([], [], color='green', marker='o', linestyle='None', markersize=8, label="Not affected by WF") 

            # Add legend
            ax.legend(handles= [gen_marker_legend] + [load_marker_legend] + [load_gen_marker_legend] + [no_load_no_gen_marker_legend] + [not_affected_legend] + [affected_legend] + edge_legend_handles, loc=loc)

        # Add colorbar for numerical values
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=np.min(self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_buses_index]), vmax=np.max(self.EF_slacks_data['contingency'][0]['node']['slack'][self.slack_buses_index])))
        sm.set_array([])  # Required for colorbar to work

        # Add colorbar that applies to both scatter plots
        cbar = plt.colorbar(sm, ax=ax, shrink=shrink, aspect=30)
        cbar.set_label('Node Slack Value')

        if title != None:
            ax.set_title(title)

        # add the California map background
        ctx.add_basemap(ax, crs=fire_area_4326_gdf.crs.to_string(), zorder=0)

        # Display the plot
        plt.show()

    # The plot showing what lines and buses were turned off
    def Off_info_plot(self, fire_area_gdf, ylim = None, xlim= None, title = None, loc='best', node_size=10):
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

        for x, y, color in zip(self.x_segments, self.y_segments, self.line_colors):
            ax.plot(x, y, color=color, linewidth=1, zorder=2)

        # Overlay the nodes with a higher zorder
        ax.scatter(self.nodes_lon, self.nodes_lat, color=self.bus_colors, s=node_size, zorder=3)

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


    def count_values_in_bins(self, vals, bin_size=10):
        """
        Manually counts the occurrences of values in bins from -100 to 100 with a given bin size.
        
        Args:
            vals (list or numpy array): Input vector of values.
            bin_size (int): The size of each bin (default is 10).
        
        Returns:
            list: Count of values in each bin.
            list: Bin labels for plotting.
        """
        # Define bin edges from -100 to 100 in increments of bin_size
        bin_edges = np.arange(-100, 100 + bin_size, bin_size)
        bin_counts = [0] * (len(bin_edges) - 1)  # Initialize count vector

        # Iterate through vals and count occurrences in each bin
        for val in vals:
            for i in range(len(bin_edges) - 1):
                if bin_edges[i] <= val < bin_edges[i + 1]:
                    bin_counts[i] += 1
                    break  # Move to next value once it's assigned to a bin

        return bin_counts, [f"{bin_edges[i]} to {bin_edges[i+1]}" for i in range(len(bin_edges)-1)]

    def plot_manual_histogram(self, vals):
        """Counts occurrences manually and plots a histogram."""
        bin_counts, bin_labels = self.count_values_in_bins(vals)
        bin_positions = np.arange(len(bin_counts)) * 10 - 100  # X-axis positions

        plt.figure(figsize=(10, 6))
        plt.bar(bin_positions, bin_counts, width=10, edgecolor='black', align='edge')
        plt.xlabel("Percentage Range")
        plt.ylabel("Frequency")
        plt.xticks(bin_positions, bin_labels, rotation=45)
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.show()

        return bin_counts
        
        print(f"Finished making con file for contingency: {Contingency_title}")

    # The power demand from every type of node
    def Power_System_demand(self):
        self.A_Power_demand = []
        self.R_Power_demand = []

        for bus in self.cats_case['bus']:
            self.A_Power_demand.append(bus['Pd'])
            self.R_Power_demand.append(bus['Qd'])

        self.A_Power_demand = np.array(self.A_Power_demand)
        self.R_Power_demand = np.array(self.R_Power_demand)

        total_Pd = sum(self.A_Power_demand)
        total_Qd = sum(self.R_Power_demand)
        l_total_Pd = sum(self.A_Power_demand[self.load_buses_only_index])
        l_total_Qd = sum(self.R_Power_demand[self.load_buses_only_index])
        g_total_Pd = sum(self.A_Power_demand[self.gen_buses_only_index])
        g_total_Qd = sum(self.R_Power_demand[self.gen_buses_only_index])
        lg_total_Pd = sum(self.A_Power_demand[self.gen_load_buses_index])
        lg_total_Qd = sum(self.R_Power_demand[self.gen_load_buses_index])

        print(f"Total System Active Power Demand: {total_Pd} MW")
        print(f"Total System Reactive Power Demand: {total_Qd} MVAR")
        print(f"Total Load only node Active Power Demand: {l_total_Pd} MW")
        print(f"Total Load only node Reactive Power Demand: {l_total_Qd} MVAR")
        print(f"Total generator only node Active Power Demand: {g_total_Pd} MW")
        print(f"Total generator only node Reactive Power Demand: {g_total_Qd} MVAR")
        print(f"Total Load and generator node Active Power Demand: {lg_total_Pd} MW")
        print(f"Total Load and generator node Reactive Power Demand: {lg_total_Qd} MVAR")

    # The power demand from every type of node that is off
    def Power_System_demand_off_nodes(self):
        if len(self.A_Power_demand ) ==0:
            self.A_Power_demand = []
            self.R_Power_demand = []

            for bus in self.cats_case['bus']:
                self.A_Power_demand.append(bus['Pd'])
                self.R_Power_demand.append(bus['Qd'])
            
            self.A_Power_demand = np.array(self.A_Power_demand)
            self.R_Power_demand = np.array(self.R_Power_demand)

        off_only_load_index = np.intersect1d(self.off_nodes_index, self.load_buses_only_index).astype(int)
        off_only_generator_index = np.intersect1d(np.array(self.off_gen_nodes)-1, self.gen_buses_only_index).astype(int)
        off_load_generator_index = np.intersect1d(np.union1d(self.off_nodes_index, np.array(self.off_gen_nodes)-1), self.gen_load_buses_index).astype(int)
        off_node_index = np.union1d(off_only_load_index, np.union1d(off_only_generator_index, off_load_generator_index)).astype(int)

        total_Pd = sum(self.A_Power_demand[off_node_index])
        total_Qd = sum(self.R_Power_demand[off_node_index])
        l_total_Pd = sum(self.A_Power_demand[off_only_load_index ])
        l_total_Qd = sum(self.R_Power_demand[off_only_load_index ])
        g_total_Pd = sum(self.A_Power_demand[off_only_generator_index ])
        g_total_Qd = sum(self.R_Power_demand[off_only_generator_index ])
        lg_total_Pd = sum(self.A_Power_demand[off_load_generator_index])
        lg_total_Qd = sum(self.R_Power_demand[off_load_generator_index])

        print(f"Total off node Active Power Demand: {total_Pd} MW")
        print(f"Total off node Reactive Power Demand: {total_Qd} MVAR")
        print(f"Total off load only node Active Power Demand: {l_total_Pd} MW")
        print(f"Total off load only node Reactive Power Demand: {l_total_Qd} MVAR")
        print(f"Total off generator only node Active Power Demand: {g_total_Pd} MW")
        print(f"Total off generator only node Reactive Power Demand: {g_total_Qd} MVAR")
        print(f"Total off load and generator node Active Power Demand: {lg_total_Pd} MW")
        print(f"Total off load and generator node Reactive Power Demand: {lg_total_Qd} MVAR")

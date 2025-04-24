
import numpy as np
import re

# Parse the data for either ACOPF solution or SCACOPF basecase solution from ExaJuGO
def parse_solution1(filename):
    with open(filename, "r") as file:
        lines = file.readlines()

    solution_data = {"bus": [], "generator": []}
    bus_data = {"i": [], "v": [], "theta": [], "bcs": []}
    generator_data = {"i": [], "id": [], "p": [], "q": []}
    
    current_section = None  # Track whether we are in bus or generator section
    skip_header = False  # Flag to skip the header row

    for line in lines:
        line = line.strip()
        
        # Identify sections
        if line == "--bus section":
            current_section = "bus"
            skip_header = True
            continue
        elif line == "--generator section":
            current_section = "generator"
            skip_header = True
            continue
        
        # Skip header row
        if skip_header:
            skip_header = False
            continue
        
        # Parse Bus Data
        if current_section == "bus":
            parts = list(map(float, line.split(",")))
            bus_data["i"].append(int(parts[0]))  # Bus index as integer
            bus_data["v"].append(parts[1])  # Voltage magnitude
            bus_data["theta"].append(parts[2])  # Angle
            bus_data["bcs"].append(parts[3])  # MVAR compensation

        # Parse Generator Data
        elif current_section == "generator":
            parts = re.split(r",\s*", line)  # Handle spaces after commas
            generator_data["i"].append(int(parts[0]))  # Bus index as integer
            generator_data["id"].append(parts[1].strip("'"))  # Remove quotes from ID
            generator_data["p"].append(float(parts[2]))  # Active power (MW)
            generator_data["q"].append(float(parts[3]))  # Reactive power (MW)

    # Convert lists to numpy arrays
    for key in bus_data:
        bus_data[key] = np.array(bus_data[key])
    for key in generator_data:
        generator_data[key] = np.array(generator_data[key], dtype=object if key == "id" else float)

    solution_data["bus"] = bus_data
    solution_data["generator"] = generator_data

    return solution_data

# Parse the data for SCACOPF Contingency solutions from ExaJuGO
def parse_solution2(filename):
    contingencies = []
    current_contingency = None
    section = None  

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line or line.startswith("--delta section") or line.startswith("delta(MW)"):  
                continue  

            if line.startswith("--contingency"):
                if current_contingency is not None:
                    contingencies.append(current_contingency)  # Save previous contingency
                current_contingency = {"label": "", "bus": {"i": [], "v": [], "theta": [], "bcs": []},
                                       "generator": {"i": [], "id": [], "p": [], "q": []}}
                section = None  
                continue

            if line.startswith("label"):
                next_line = next(file, "").strip()  
                if next_line:  
                    current_contingency["label"] = next_line.strip("'")  
                continue

            if line.startswith("--bus section"):
                section = "bus"
                continue

            if line.startswith("--generator section"):
                section = "generator"
                continue

            if section == "bus" and line and not line.startswith("i"):
                parts = line.split(", ")
                if len(parts) == 4:
                    current_contingency["bus"]["i"].append(int(parts[0]))
                    current_contingency["bus"]["v"].append(float(parts[1]))
                    current_contingency["bus"]["theta"].append(float(parts[2]))
                    current_contingency["bus"]["bcs"].append(float(parts[3]))

            elif section == "generator" and line and not line.startswith("i"):
                parts = line.split(", ")
                if len(parts) == 4:
                    current_contingency["generator"]["i"].append(int(parts[0]))
                    current_contingency["generator"]["id"].append(parts[1].strip("'"))  
                    current_contingency["generator"]["p"].append(float(parts[2]))
                    current_contingency["generator"]["q"].append(float(parts[3]))

        if current_contingency is not None:
            contingencies.append(current_contingency)  

    return contingencies

# Parse the power flow data from Exajugo ACOPF or SCACOPF
def parse_power_data(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    
    base_data = {"from": [], "to": [], "p_from": [], "p_to": []}
    contingency_cases = []  # List of independent contingency cases
    current_contingency = None
    section = None  # Keeps track of whether we are in transformer/line sections
    is_contingency = False  # Flag for contingency sections
    current_label = None  # Store label of the contingency

    for line in lines:
        line = line.strip()
        
        if line == "--base":
            is_contingency = False
        elif line == "--contingency":
            is_contingency = True
            section = None  # Reset section indicator
        elif line == "label":
            continue
        elif line.startswith("'"):
            # Start a new contingency case
            current_label = line.strip("'")  # Save the label
            if current_contingency is not None:
                contingency_cases.append(current_contingency)  # Save previous case
            current_contingency = {"label": current_label, "from": [], "to": [], "p_from": [], "p_to": []}
        elif "transformer section" in line or "line section" in line:
            section = "data"  # Switch to data mode
        elif section == "data" and line and not line.startswith("From"):
            # Extract numeric data
            parts = line.split(", ")
            f, t, p_f, p_t = int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3])
            
            target = current_contingency if is_contingency else base_data
            target["from"].append(f)
            target["to"].append(t)
            target["p_from"].append(p_f)
            target["p_to"].append(p_t)
    
    # Append the last contingency case if it exists
    if current_contingency is not None:
        contingency_cases.append(current_contingency)

    # Convert lists to numpy arrays
    for key in base_data:
        base_data[key] = np.array(base_data[key])
    for case in contingency_cases:
        for key in ["from", "to", "p_from", "p_to"]:
            case[key] = np.array(case[key])

    return base_data, contingency_cases   

# Parse the components of the power flow constraints from Exajugo ACOPF or SCACOPF
def parse_power_constraint_data(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    
    p_data = {"Bus_id": [], "sum_pg": [], "Pd": [], "Gsh": [], "v_n": [], "Gvn2": [], "sum_p_li": [], 
                "sum_p_ti": [], "pslackp_n": [], "pslackm_n": [], "p_relax_slack": []}
    q_data = {"Bus_id": [], "sum_qg": [], "Qd": [], "Bsh": [], "b_s": [], "v_n": [], "Bmb_svn2": [], 
                "sum_q_li": [], "sum_q_ti": [], "qslackp_n": [], "qslackm_n": [], "q_relax_slack": []}
    
    current_section = None  # Track the current section
    skip_header = False  # Flag to skip header row
    
    for line in lines:
        line = line.strip()
        
        # Identify sections
        if line == "--Active power constraint":
            current_section = "Active"
            skip_header = True
            continue
        elif line == "--Reactive power constraint":
            current_section = "Reactive"
            skip_header = True
            continue
        
        # Skip header row
        if skip_header:
            skip_header = False
            continue
        
        # Parse Node Slack Data
        if current_section == "Active":
            parts = list(map(float, line.split(",")))
            p_data["Bus_id"].append(int(parts[0]))
            p_data["sum_pg"].append(parts[1])
            p_data["Pd"].append(parts[2])
            p_data["Gsh"].append(parts[3])
            p_data["v_n"].append(parts[4])
            p_data["Gvn2"].append(parts[5])
            p_data["sum_p_li"].append(parts[6])
            p_data["sum_p_ti"].append(parts[7])
            p_data["pslackm_n"].append(parts[8])
            p_data["pslackp_n"].append(parts[9])
            p_data["p_relax_slack"].append(parts[10])
        
        # Parse Line Slack Data
        elif current_section == "Reactive":
            parts = list(map(float, line.split(",")))
            q_data["Bus_id"].append(int(parts[0]))
            q_data["sum_qg"].append(parts[1])
            q_data["Qd"].append(parts[2])
            q_data["Bsh"].append(parts[3])
            q_data["b_s"].append(parts[4])
            q_data["v_n"].append(parts[5])
            q_data["Bmb_svn2"].append(parts[6])
            q_data["sum_q_li"].append(parts[7])
            q_data["sum_q_ti"].append(parts[8])
            q_data["qslackm_n"].append(parts[9])
            q_data["qslackp_n"].append(parts[10])
            q_data["q_relax_slack"].append(parts[11])
    
    # Convert lists to numpy arrays
    for key in p_data:
        p_data[key] = np.array(p_data[key])
    for key in q_data:
        q_data[key] = np.array(q_data[key])
    
    return p_data, q_data

# Parse the slack data for the buses and lines from Exajugo ACOPF or SCACOPF basecase
def parse_slack_data(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    
    basecase = {"node_slack": {}, "line_slack": {}, "transformer_slack": {}}
    node_slack_data = {"i": [], "pslackm_n": [], "pslackp_n": [], "qslackm_n": [], "qslackp_n": []}
    line_slack_data = {"from": [], "to": [], "sslack_li_from": [], "sslack_li_to": []}
    transformer_slack_data = {"from": [], "to": [], "sslack_ti_from": [], "sslack_ti_to": []}
    
    current_section = None  # Track the current section
    skip_header = False  # Flag to skip header row
    
    for line in lines:
        line = line.strip()
        
        # Identify sections
        if line == "--node slack section":
            current_section = "node"
            skip_header = True
            continue
        elif line == "--line slack section":
            current_section = "line"
            skip_header = True
            continue
        elif line == "--transformer slack section":
            current_section = "transformer"
            skip_header = True
            continue
        
        # Skip header row
        if skip_header:
            skip_header = False
            continue
        
        # Parse Node Slack Data
        if current_section == "node":
            parts = list(map(float, line.split(",")))
            node_slack_data["i"].append(int(parts[0]))
            node_slack_data["pslackm_n"].append(parts[1])
            node_slack_data["pslackp_n"].append(parts[2])
            node_slack_data["qslackm_n"].append(parts[3])
            node_slack_data["qslackp_n"].append(parts[4])
        
        # Parse Line Slack Data
        elif current_section == "line":
            parts = list(map(float, line.split(",")))
            line_slack_data["from"].append(int(parts[0]))
            line_slack_data["to"].append(int(parts[1]))
            line_slack_data["sslack_li_from"].append(parts[2])
            line_slack_data["sslack_li_to"].append(parts[3])
        
        # Parse Transformer Slack Data
        elif current_section == "transformer":
            parts = list(map(float, line.split(",")))
            transformer_slack_data["from"].append(int(parts[0]))
            transformer_slack_data["to"].append(int(parts[1]))
            transformer_slack_data["sslack_ti_from"].append(parts[2])
            transformer_slack_data["sslack_ti_to"].append(parts[3])
    
    # Convert lists to numpy arrays
    for key in node_slack_data:
        node_slack_data[key] = np.array(node_slack_data[key])
    for key in line_slack_data:
        line_slack_data[key] = np.array(line_slack_data[key])
    for key in transformer_slack_data:
        transformer_slack_data[key] = np.array(transformer_slack_data[key])
    
    basecase["node_slack"] = node_slack_data
    basecase["line_slack"] = line_slack_data
    basecase["transformer_slack"] = transformer_slack_data
    return basecase

# Parse the slack data for the buses and lines from Exajugo SCACOPF contingencies
def parse_SCACOPF_slack(filename):
    contingencies = []
    current_contingency = None
    section = None  

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith("--contingency"):
                if current_contingency is not None:
                    contingencies.append(current_contingency)  # Save previous contingency
                current_contingency = {"label": "", 
                                       "node_slack": {"i": [], "pslackm_n": [], "pslackp_n": [], "qslackm_n": [], "qslackp_n": []},
                                       "line_slack": {"from": [], "to": [], "sslack_li_from": [], "sslack_li_to": []},
                                       "transformer_slack": {"from": [], "to": [], "sslack_ti_from": [], "sslack_ti_to": []}}
                section = None  
                continue
            elif line.startswith("--base"):
                current_contingency = {"label": "Basecase", 
                                       "node_slack": {"i": [], "pslackm_n": [], "pslackp_n": [], "qslackm_n": [], "qslackp_n": []},
                                       "line_slack": {"from": [], "to": [], "sslack_li_from": [], "sslack_li_to": []},
                                       "transformer_slack": {"from": [], "to": [], "sslack_ti_from": [], "sslack_ti_to": []}}
                section = None  
                continue

            if line.startswith("--label"):
                next_line = next(file, "").strip() 
                if next_line:  
                    current_contingency["label"] = next_line.strip("'")  
                continue

            if line.startswith("--node slack section"):
                section = "node_slack"
                continue

            if line.startswith("--line slack section"):
                section = "line_slack"
                continue

            if line.startswith("--transformer slack section"):
                section = "transformer_slack"
                continue

            if section == "node_slack" and line and not line.startswith("Bus"):
                parts = line.split(", ")
                if len(parts) == 5:
                    current_contingency["node_slack"]["i"].append(int(parts[0]))
                    current_contingency["node_slack"]["pslackm_n"].append(float(parts[1]))
                    current_contingency["node_slack"]["pslackp_n"].append(float(parts[2]))
                    current_contingency["node_slack"]["qslackm_n"].append(float(parts[3]))
                    current_contingency["node_slack"]["qslackp_n"].append(float(parts[4]))

            elif section == "line_slack" and line and not line.startswith("From"):
                parts = line.split(", ")
                if len(parts) == 4:
                    current_contingency["line_slack"]["from"].append(int(parts[0]))
                    current_contingency["line_slack"]["to"].append(int(parts[1]))
                    current_contingency["line_slack"]["sslack_li_from"].append(float(parts[2]))
                    current_contingency["line_slack"]["sslack_li_to"].append(float(parts[3]))

            elif section == "transformer_slack" and line and not line.startswith("From"):
                parts = line.split(", ")
                if len(parts) == 4:
                    current_contingency["transformer_slack"]["from"].append(int(parts[0]))
                    current_contingency["transformer_slack"]["to"].append(int(parts[1]))
                    current_contingency["transformer_slack"]["sslack_ti_from"].append(float(parts[2]))
                    current_contingency["transformer_slack"]["sslack_ti_to"].append(float(parts[3]))

        if current_contingency is not None:
            contingencies.append(current_contingency)  

        for i in range(len(contingencies)):
            for key in contingencies[i]['node_slack']:
                contingencies[i]['node_slack'][key] = np.array(contingencies[i]['node_slack'][key])
            for key in contingencies[i]['line_slack']:
                contingencies[i]['line_slack'][key] = np.array(contingencies[i]['line_slack'][key])
            for key in contingencies[i]['transformer_slack']:
                contingencies[i]['transformer_slack'][key] = np.array(contingencies[i]['transformer_slack'][key])

    return contingencies

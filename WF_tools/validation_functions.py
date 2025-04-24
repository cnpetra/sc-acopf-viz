from parse_data import *

# Check that the slack for the lines and buses below than a negative bound
def check_slacks(slacks, neg_bound):
    pmn_large_neg = np.where(slacks['node_slack']['pslackm_n'] < neg_bound)[0]
    if len(pmn_large_neg) > 0:
        print('Bus ' + str(slacks['node_slack']['i'][pmn_large_neg]) 
                + ' have large negative slack for pslackm_n')

    ppn_large_neg = np.where(slacks['node_slack']['pslackp_n'] < neg_bound)[0]
    if len(ppn_large_neg) > 0:
        print('Bus ' + str(slacks['node_slack']['i'][ppn_large_neg]) 
                + ' have large negative slack for pslackp_n')

    qmn_large_neg = np.where(slacks['node_slack']['qslackm_n'] < neg_bound)[0]
    if len(qmn_large_neg) > 0:
        print('Bus ' + str(slacks['node_slack']['i'][qmn_large_neg]) 
                + ' have large negative slack for qslackm_n')

    qpn_large_neg = np.where(slacks['node_slack']['qslackp_n'] < neg_bound)[0]
    if len(qpn_large_neg) > 0:
        print('Bus ' + str(slacks['node_slack']['i'][qpn_large_neg]) 
                + ' have large negative slack for qslackp_n')

    sfrom_large_neg = np.where(slacks['line_slack']['sslack_li_from'] < neg_bound)[0]
    if len(sfrom_large_neg) > 0:
        print('Line with from bus and to bus of ' + str(slacks['line_slack']['from'][sfrom_large_neg]) 
                + ' '  + str(slacks['line_slack']['to'][sfrom_large_neg]) 
                + ' have large negative slack for sslack_li_from')

    sto_large_neg = np.where(slacks['line_slack']['sslack_li_to'] < neg_bound)[0]
    if len(sto_large_neg) > 0:
        print('Line with from bus and to bus of ' + str(slacks['line_slack']['from'][sto_large_neg]) 
                + ' '  + str(slacks['line_slack']['to'][sto_large_neg]) 
                + ' have large negative slack for sslack_li_to')

    sfrom_large_neg = np.where(slacks['transformer_slack']['sslack_ti_from'] < neg_bound)[0]
    if len(sfrom_large_neg) > 0:
        print('Line with from bus and to bus of ' + str(slacks['transformer_slack']['from'][sfrom_large_neg]) 
                + ' '  + str(slacks['transformer_slack']['to'][sfrom_large_neg]) 
                + ' have large negative slack for sslack_ti_from')

    sto_large_neg = np.where(slacks['transformer_slack']['sslack_ti_to'] < neg_bound)[0]
    if len(sto_large_neg) > 0:
        print(len(sto_large_neg))
        print(sto_large_neg)
        print('Line with from bus and to bus of ' + str(slacks['transformer_slack']['from'][sto_large_neg]) 
                + ' '  + str(slacks['transformer_slack']['to'][sto_large_neg]) 
                + ' have large negative slack for sslack_ti_to')

# Collect the slack data for the iterative from of SCACOPF problem 
# The iterative form of SCACOPF solves the basecase and contingency separely  
def collect_IF_slack_data(ACOPF_slack_filename, Contingency_slack_filename):

    IF_bc_slacks = parse_slack_data(ACOPF_slack_filename)
    IF_con_slacks = parse_slack_data(Contingency_slack_filename)
    IF_slack_data = {'basecase': [], 'contingency': []}
    basecase = {'node':[], 'line': [], 'transformer': []}
    contingencies = []
    contingency_k = {'node':[], 'line': [], 'transformer': []}

    # Start the basecase of the data
    node_data = {'i': IF_bc_slacks['node_slack']['i'], 'slack': []}
    line_data = {'from': IF_bc_slacks['line_slack']['from'], 'to': IF_bc_slacks['line_slack']['to'], 'slack': []}
    transformer_data = {'from': IF_bc_slacks['transformer_slack']['from'], 'to': IF_bc_slacks['transformer_slack']['to'], 'slack': []}

    IF_bc_bus_slack = np.maximum.reduce([IF_bc_slacks['node_slack']['pslackm_n'], IF_bc_slacks['node_slack']['pslackp_n'], 
                                    IF_bc_slacks['node_slack']['qslackm_n'], IF_bc_slacks['node_slack']['qslackp_n']])

    IF_bc_line_slack = np.maximum.reduce([IF_bc_slacks['line_slack']['sslack_li_from'], 
                                    IF_bc_slacks['line_slack']['sslack_li_to']])

    IF_bc_transformer_slack = np.maximum.reduce([IF_bc_slacks['transformer_slack']['sslack_ti_from'], 
                                            IF_bc_slacks['transformer_slack']['sslack_ti_to']])

    check_slacks(IF_bc_slacks, -1)

    IF_bc_bus_slack[IF_bc_bus_slack < 0] = 0
    IF_bc_line_slack[IF_bc_line_slack < 0] = 0
    IF_bc_transformer_slack[IF_bc_transformer_slack < 0] = 0

    node_data['slack'] = IF_bc_bus_slack
    line_data['slack'] = IF_bc_line_slack
    transformer_data['slack'] = IF_bc_transformer_slack

    basecase['node'] = node_data
    basecase['line'] = line_data
    basecase['transformer'] = transformer_data

    # Start the contingency part of the data
    node_data = {'i': IF_con_slacks['node_slack']['i'], 'slack': []}
    line_data = {'from': IF_con_slacks['line_slack']['from'], 'to': IF_con_slacks['line_slack']['to'], 'slack': []}
    transformer_data = {'from': IF_con_slacks['transformer_slack']['from'], 'to': IF_con_slacks['transformer_slack']['to'], 'slack': []}

    IF_con_bus_slack = np.maximum.reduce([IF_con_slacks['node_slack']['pslackm_n'], IF_con_slacks['node_slack']['pslackp_n'], 
                                    IF_con_slacks['node_slack']['qslackm_n'], IF_con_slacks['node_slack']['qslackp_n']])

    IF_con_line_slack = np.maximum.reduce([IF_con_slacks['line_slack']['sslack_li_from'], 
                                    IF_con_slacks['line_slack']['sslack_li_to']])

    IF_con_transformer_slack = np.maximum.reduce([IF_con_slacks['transformer_slack']['sslack_ti_from'], 
                                            IF_con_slacks['transformer_slack']['sslack_ti_to']])

    check_slacks(IF_con_slacks, -1)

    IF_con_bus_slack[IF_con_bus_slack < 0] = 0
    IF_con_line_slack[IF_con_line_slack < 0] = 0
    IF_con_transformer_slack[IF_con_transformer_slack < 0] = 0

    node_data['slack'] = IF_con_bus_slack
    line_data['slack'] = IF_con_line_slack
    transformer_data['slack'] = IF_con_transformer_slack

    contingency_k['node'] = node_data
    contingency_k['line'] = line_data
    contingency_k['transformer'] = transformer_data

    contingencies.append(contingency_k)

    # put everything together
    IF_slack_data['basecase'] = basecase
    IF_slack_data['contingency'] = contingencies

    return IF_slack_data

# Collect the slack data for the extended from of SCACOPF problem 
# The extended form of SCACOPF solves the basecase and contingency as one big optimization problem  
def collect_EF_slack_data(SCACOPF_slack_filename):

    EF_slacks = parse_SCACOPF_slack(SCACOPF_slack_filename)

    EF_slack_data = {'basecase': [], 'contingency': []}
    basecase = {'node':[], 'line': [], 'transformer': []}
    contingencies = []
    contingency_k = {'node':[], 'line': [], 'transformer': []}

    # Start the basecase of the data
    node_data = {'i': EF_slacks[0]['node_slack']['i'], 'slack': []}
    line_data = {'from': EF_slacks[0]['line_slack']['from'], 'to': EF_slacks[0]['line_slack']['to'], 'slack': []}
    transformer_data = {'from': EF_slacks[0]['transformer_slack']['from'], 'to': EF_slacks[0]['transformer_slack']['to'], 'slack': []}

    EF_bc_bus_slack = np.maximum.reduce([EF_slacks[0]['node_slack']['pslackm_n'], EF_slacks[0]['node_slack']['pslackp_n'], 
                                        EF_slacks[0]['node_slack']['qslackm_n'], EF_slacks[0]['node_slack']['qslackp_n']])

    EF_bc_line_slack = np.maximum.reduce([EF_slacks[0]['line_slack']['sslack_li_from'], 
                                            EF_slacks[0]['line_slack']['sslack_li_to']])

    EF_bc_transformer_slack = np.maximum.reduce([EF_slacks[0]['transformer_slack']['sslack_ti_from'], 
                                                    EF_slacks[0]['transformer_slack']['sslack_ti_to']])

    check_slacks(EF_slacks[0], -1)

    EF_bc_bus_slack[EF_bc_bus_slack < 0] = 0
    EF_bc_line_slack[EF_bc_line_slack < 0] = 0
    EF_bc_transformer_slack[EF_bc_transformer_slack < 0] = 0

    node_data['slack'] = EF_bc_bus_slack
    line_data['slack'] = EF_bc_line_slack
    transformer_data['slack'] = EF_bc_transformer_slack

    basecase['node'] = node_data
    basecase['line'] = line_data
    basecase['transformer'] = transformer_data

    # Start the contingency part of the data
    node_data = {'i': EF_slacks[1]['node_slack']['i'], 'slack': []}
    line_data = {'from': EF_slacks[1]['line_slack']['from'], 'to': EF_slacks[1]['line_slack']['to'], 'slack': []}
    transformer_data = {'from': EF_slacks[1]['transformer_slack']['from'], 'to': EF_slacks[1]['transformer_slack']['to'], 'slack': []}

    EF_con_bus_slack = np.maximum.reduce([EF_slacks[1]['node_slack']['pslackm_n'], EF_slacks[1]['node_slack']['pslackp_n'], 
                                        EF_slacks[1]['node_slack']['qslackm_n'], EF_slacks[1]['node_slack']['qslackp_n']])

    EF_con_line_slack = np.maximum.reduce([EF_slacks[1]['line_slack']['sslack_li_from'], 
                                            EF_slacks[1]['line_slack']['sslack_li_to']])

    EF_con_transformer_slack = np.maximum.reduce([EF_slacks[1]['transformer_slack']['sslack_ti_from'], 
                                                    EF_slacks[1]['transformer_slack']['sslack_ti_to']])

    check_slacks(EF_slacks[1], -1)

    EF_con_bus_slack[EF_con_bus_slack < 0] = 0
    EF_con_line_slack[EF_con_line_slack < 0] = 0
    EF_con_transformer_slack[EF_con_transformer_slack < 0] = 0

    node_data['slack'] = EF_con_bus_slack 
    line_data['slack'] = EF_con_line_slack
    transformer_data['slack'] = EF_con_transformer_slack

    contingency_k['node'] = node_data
    contingency_k['line'] = line_data
    contingency_k['transformer'] = transformer_data

    contingencies.append(contingency_k)

    # put everything together
    EF_slack_data['basecase'] = basecase
    EF_slack_data['contingency'] = contingencies

    return EF_slack_data
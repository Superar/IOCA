#!/usr/bin/python
# -*- coding: utf-8 -*-

from collections import namedtuple
from ortools.linear_solver import pywraplp
import math

Point = namedtuple("Point", ['x', 'y'])
Facility = namedtuple("Facility", ['index', 'setup_cost', 'capacity', 'location'])
Customer = namedtuple("Customer", ['index', 'demand', 'location'])

DEBUG = 0

def length(point1, point2):
    return math.sqrt((point1.x - point2.x)**2 + (point1.y - point2.y)**2)

def solve_it(input_data):
    # parse the input
    lines = input_data.split('\n')

    parts = lines[0].split()
    facility_count = int(parts[0])
    customer_count = int(parts[1])

    facilities = []
    for i in range(1, facility_count+1):
        parts = lines[i].split()
        facilities.append(Facility(i-1, float(parts[0]), int(parts[1]), Point(float(parts[2]), float(parts[3])) ))

    customers = []
    for i in range(facility_count+1, facility_count+1+customer_count):
        parts = lines[i].split()
        customers.append(Customer(i-1-facility_count, int(parts[0]), Point(float(parts[1]), float(parts[2]))))

    return facilityNaive(facility_count, facilities, customer_count, customers)


def facilityNaive(facility_count, facilities, customer_count, customers):

    if DEBUG >= 1:
        print(f"Numero de possiveis instalacoes = {facility_count}")
        print(f"Numero de clientes = {customer_count}")

    if DEBUG >= 2:
        print("Instalacoes na ordem que foram lidas")
        for facility in facilities:
            print(facility)
        print()

    if DEBUG >= 2:
        print("Clientes na ordem que foram lidos")
        for customer in customers:
            print(customer)
        print()

    # Modify this code to run your optimization algorithm
    solutions = list()

    # -------------------------------------------------- #
    # trivial solution: pack the facilities one by one until all the customers are served
    solution_trivial = [-1]*len(customers)
    capacity_remaining = [f.capacity for f in facilities]
    facility_index = 0
    for customer in customers:
        if capacity_remaining[facility_index] >= customer.demand:
            solution_trivial[customer.index] = facility_index
            capacity_remaining[facility_index] -= customer.demand
        else:
            facility_index += 1
            assert capacity_remaining[facility_index] >= customer.demand
            solution_trivial[customer.index] = facility_index
            capacity_remaining[facility_index] -= customer.demand
    solutions.append(solution_trivial)
    # -------------------------------------------------- #

    # ILP
    solver = pywraplp.Solver('FL',
                             pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    solver.set_time_limit(120000) # 2 minute time limit
    open_vars = [solver.IntVar(0, 1, 'f[{}]'.format(f.index))
                 for f in facilities]
    connection_vars = [[solver.IntVar(0, 1, 'c[{}]--f[{}]'.format(c.index, f.index))
                      for f in facilities] for c in customers]

    for customer in customers:
        # Each customer is connected to exactly one facility
        connection_ct = solver.Constraint(1, 1, 'c_ct[{}]'.format(customer.index))
        for facility in facilities:
            connection_ct.SetCoefficient(connection_vars[customer.index][facility.index], 1)

            # Customers can only be connected to open facilities
            open_ct = solver.Constraint(0, 1, 'op_ct[{}--{}]'.format(customer.index, facility.index))
            open_ct.SetCoefficient(open_vars[facility.index], 1)
            open_ct.SetCoefficient(connection_vars[customer.index][facility.index], -1)

    # Mantain demand under facilities capacity
    for facility in facilities:
        capacity_ct = solver.Constraint(0, facility.capacity, 'cap_ct[{}]'.format(facility.index))
        for customer in customers:
            capacity_ct.SetCoefficient(connection_vars[customer.index][facility.index], customer.demand)

    objective = solver.Objective()
    objective.SetMinimization()
    # Connecting distance
    for customer in customers:
        for facility in facilities:
            dist = length(customer.location, facility.location)
            objective.SetCoefficient(connection_vars[customer.index][facility.index], dist)
    
    # Opening cost
    for facility in facilities:
        objective.SetCoefficient(open_vars[facility.index], facility.setup_cost)
    status = solver.Solve()

    if status == solver.FEASIBLE or status == solver.OPTIMAL:
        solution_ilp = [[int(f.solution_value()) for f in c].index(1) for c in connection_vars]
        solutions.append(solution_ilp)
    # -------------------------------------------------- #

    # Greedy by customer
    sorted_customers = sorted(customers,
                              key=lambda x: x.demand,
                              reverse=True)
    solution_greedy_customer = [-1]*len(customers)
    cur_capacities = [f.capacity for f in facilities]
    facility_customer_dists = [[length(c.location, f.location) for f in facilities] for c in customers]
    for customer in sorted_customers:
        dist = facility_customer_dists[customer.index]
        used = {f for f in solution_greedy_customer if f >= 0}
        # Set costumer to closest used facility
        if used:
            best_facility = None
            cur_obj = float('inf')
            for f in used:
                if dist[f] < cur_obj:
                    # Facility is better and no setup cost is needed (because it is already being used)
                    if cur_capacities[f] >= customer.demand:
                        # There is enough room
                        best_facility = f
                        cur_obj = dist[f]
            if best_facility is not None:
                solution_greedy_customer[customer.index] = best_facility
                cur_capacities[best_facility] -= customer.demand
                continue

        # Set customer to closest unused facility
        best_facility = None
        cur_obj = float('inf')

        # Find best facility that minimizes the objective
        for f in facilities:
            if f.setup_cost + dist[f.index] < cur_obj and cur_capacities[f.index] >= customer.demand:
                    best_facility = f.index
                    cur_obj = f.setup_cost + dist[f.index]
        solution_greedy_customer[customer.index] = best_facility
        cur_capacities[best_facility] -= customer.demand
    solutions.append(solution_greedy_customer)
    # -------------------------------------------------- #

    # Greedy by facility
    solution_greedy_facility = [-1]*len(customers)
    
    while -1 in solution_greedy_facility:
        unused_facilities = [f for f in facilities if f.index not in solution_greedy_facility]

        # Get which customers each facility would attract if it were open
        facility_attraction = list()
        for f in unused_facilities:
            closest_customers = sorted(customers, key=lambda c: facility_customer_dists[c.index][f.index])
            cur_capacity = f.capacity
            cur_f_attraction = list()
            for c in closest_customers:
                if solution_greedy_facility[c.index] != -1:
                    # Customer is already beeing attended
                    # See if it is worth moving it
                    cur_c_facility = facilities[solution_greedy_facility[c.index]]
                    cur_c_cost = cur_c_facility.setup_cost + facility_customer_dists[c.index][cur_c_facility.index]
                    new_c_cost = f.setup_cost + facility_customer_dists[c.index][f.index]
                    if new_c_cost < cur_c_cost and c.demand <= cur_capacity:
                        cur_capacity -= c.demand
                        cur_f_attraction.append(c)
                elif c.demand <= cur_capacity:
                    cur_capacity -= c.demand
                    cur_f_attraction.append(c)
            facility_attraction.append(cur_f_attraction)

        facility_obj = [f.setup_cost + sum([facility_customer_dists[c.index][f.index]
                                            for c in facility_attraction[i]])
                        for i, f in enumerate(unused_facilities)]
        best_facility_idx = facility_obj.index(min(facility_obj))
        for c in facility_attraction[best_facility_idx]:
            solution_greedy_facility[c.index] = unused_facilities[best_facility_idx].index
    solutions.append(solution_greedy_facility)
    # -------------------------------------------------- #

    # Local heuristic search
    initial_solution = solution_greedy_facility.copy()

    while True:
        used_facilities = {f for f in facilities if f.index in initial_solution}
        unused_facilities = {f for f in facilities if f.index not in initial_solution}
        
        facilities_capacities = [f.capacity - sum([customers[c].demand
                                                for c, f_idx in enumerate(initial_solution)
                                                if f_idx == f.index])
                                for f in facilities]
        facility_gains = dict()
        facility_attraction = dict()
        for f in unused_facilities:
            closest_customers = sorted(customers, key=lambda c: facility_customer_dists[c.index][f.index])
            cur_capacity = facilities_capacities[f.index]
            cur_f_attraction = list()
            for c in closest_customers:
                cur_c_facility = facilities[initial_solution[c.index]]
                cur_c_cost = facility_customer_dists[c.index][cur_c_facility.index]
                new_c_cost = facility_customer_dists[c.index][f.index]
                if new_c_cost < cur_c_cost and c.demand <= cur_capacity:
                    cur_capacity -= c.demand
                    cur_f_attraction.append(c)
            facility_attraction[f.index] = cur_f_attraction
            
            cur_f_gain = sum([facility_customer_dists[c.index][initial_solution[c.index]]
                            for c in cur_f_attraction])
            cur_f_gain -= sum([facility_customer_dists[c.index][f.index] for c in cur_f_attraction])
            facility_gains[f.index] = cur_f_gain
        
        customers_should_go = dict()
        for f in used_facilities:
            cur_f_gain = 0
            cur_capacity = facilities_capacities.copy()
            for c in customers:
                if initial_solution[c.index] == f.index:
                    # If facility f werent used, where would customers go?
                    cur_c_dists = facility_customer_dists[c.index]
                    closest_facilities = sorted(cur_c_dists)
                    new_position = f.index # If no other facility can receive the customer, it should stay where it is
                    for min_dist in closest_facilities:
                        if cur_c_dists.index(min_dist) != f.index and c.demand <= cur_capacity[cur_c_dists.index(min_dist)]:
                            new_position = cur_c_dists.index(min_dist)
                            cur_capacity[new_position] -= c.demand
                            break
                    cur_f_gain += facility_customer_dists[c.index][new_position]
                    cur_f_gain -= facility_customer_dists[c.index][f.index]

                    customers_should_go[c.index] = new_position
            facility_gains[f.index] = cur_f_gain

        should_use = {f.index for f in unused_facilities if f.setup_cost < facility_gains[f.index]}
        # If the gain of a used facility is 0, it means that there was no other facility where its customers could go
        # It should stay as it is
        should_not_use = {f.index for f in used_facilities if f.setup_cost > facility_gains[f.index] and facility_gains[f.index] > 0}
        change_state = should_use | should_not_use

        if change_state:
            changing_facility_idx = change_state.pop()
            print(changing_facility_idx)
            changing_facility = facilities[changing_facility_idx]
            if changing_facility_idx in should_use:
                # Set facility to those customers, which it attracts
                for c in facility_attraction[changing_facility.index]:
                    initial_solution[c.index] = changing_facility.index
            else:
                for i, f in enumerate(initial_solution):
                    if f == changing_facility.index:
                        initial_solution[i] = customers_should_go[i]
        else:
            break
    solutions.append(initial_solution)
    # -------------------------------------------------- #

    def cost(solution):
        used = [0]*len(facilities)
        for facility_index in solution:
            used[facility_index] = 1

        # calculate the cost of the solution
        obj = sum([f.setup_cost*used[f.index] for f in facilities])
        for customer in customers:
            obj += length(customer.location, facilities[solution[customer.index]].location)
        
        return obj

    solution_costs = [cost(sol) for sol in solutions]
    print(solution_costs)
    best_idx = solution_costs.index(min(solution_costs))
    obj = solution_costs[best_idx]
    solution = solutions[best_idx]

    # prepare the solution in the specified output format
    output_data = '%.2f' % obj + '\n'
    output_data += ' '.join(map(str, solution))

    return output_data


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        output_data = solve_it(input_data)
        print(output_data)
        solution_file = open(file_location + ".sol", "w")
        solution_file.write(output_data)
        solution_file.close()
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/fl_16_2)')

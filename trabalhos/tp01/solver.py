from collections import namedtuple
from ortools.linear_solver import pywraplp

Item = namedtuple("Item", ['index', 'value', 'weight'])

DEBUG = 0


def solve_it(input_data):
    # parse the input
    lines = input_data.split('\n')

    firstLine = lines[0].split()
    item_count = int(firstLine[0])
    capacity = int(firstLine[1])
    conflict_count = int(firstLine[2])

    items = []
    conflicts = []

    for i in range(1, item_count+1):
        line = lines[i]
        parts = line.split()
        items.append(Item(i-1, int(parts[0]), int(parts[1])))

    for i in range(1, conflict_count+1):
        line = lines[item_count + i]
        parts = line.split()
        conflicts.append((int(parts[0]), int(parts[1])))

    return knapsackNaive(item_count, items, capacity, conflict_count, conflicts)


def get_neighbours(cur_solution, items, capacity, conflicts):
    neighbours = list()
    for i, value in enumerate(cur_solution):
        neighbour = cur_solution.copy()
        neighbour[i] += pow(-1, value)
        neighbours.append(neighbour)

    # Remove solutions which exceed capacity
    neighbours_weights = [sum([items[i].weight for i, v in enumerate(n) if v])
                          for n in neighbours]
    neighbours = [n for i, n in enumerate(neighbours)
                  if neighbours_weights[i] <= capacity]

    # Remove solutions with conflicts
    for i, j in conflicts:
        conflict_neighbours = [k for k, n in enumerate(neighbours)
                               if n[i] and n[j]]
        neighbours = [n for k, n in enumerate(neighbours)
                      if k not in conflict_neighbours]
    return neighbours


def knapsackNaive(num_items, items, capacity, num_conflicts, conflicts):

    if DEBUG >= 1:
        print(f"numero de itens = {num_items}")
        print(f"capacidade da mochila = {capacity}")
        print(f"numero de conflitos = {num_conflicts}")

    if DEBUG >= 2:
        print("Itens na ordem em que foram lidos")
        for item in items:
            print(item)
        print()

    if DEBUG >= 2:
        print("Conflitos na ordem em que foram lidos")
        for conflict in conflicts:
            print(conflict)
        print()

    # Modify this code to run your optimization algorithm

    # ILP
    solver = pywraplp.Solver('MochilaColorida',
                             pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
    pli_vars = [solver.IntVar(0, 1, 'x[{}]'.format(i))
                for i, _ in enumerate(items)]
    ct = solver.Constraint(0, capacity, 'capacity')
    for i, item in enumerate(items):
        ct.SetCoefficient(pli_vars[i], item.weight)

    for i, (item1, item2) in enumerate(conflicts):
        item_ct = solver.Constraint(0, 1, 'item_ct[{}]'.format(i))
        item_ct.SetCoefficient(pli_vars[item1], 1)
        item_ct.SetCoefficient(pli_vars[item2], 1)

    objective = solver.Objective()
    objective.SetMaximization()
    for i, item in enumerate(items):
        objective.SetCoefficient(pli_vars[i], item.value)
    solver.Solve()

    solution_value_ilp = int(objective.Value())
    solution_ilp = [int(x.solution_value()) for x in pli_vars]

    # Greedy
    solution_greedy = [0]*num_items
    solution_value_greedy = 0
    cur_weight = 0
    custo_beneficio = [item.value*1.0/item.weight for item in items]
    sorted_items = sorted(items, key=lambda x: custo_beneficio[x.index])
    banned_items = list()

    while cur_weight <= capacity and sorted_items:
        item = sorted_items.pop()
        if cur_weight + item.weight <= capacity and item.index not in banned_items:
            cur_weight += item.weight
            solution_greedy[item.index] = 1
            solution_value_greedy += item.value
            if num_conflicts:
                banned_items.extend(
                    [i2 for i1, i2 in conflicts if i1 == item.index])
                banned_items.extend(
                    [i1 for i1, i2 in conflicts if i2 == item.index])
        else:
            continue

    # Local search heuristic
    old_value = -1
    solution_search = solution_greedy
    solution_value_search = solution_value_greedy

    while solution_value_search > old_value:
        neighbours = get_neighbours(
            solution_search, items, capacity, conflicts)
        neighbours_values = [sum([items[i].value for i, v in enumerate(n) if v])
                             for n in neighbours]
        max_neighbour_value = max(neighbours_values)

        if max_neighbour_value > solution_value_search:
            old_value = solution_value_search
            solution_value_search = max_neighbour_value
            max_neighbour_index = neighbours_values.index(max_neighbour_value)
            solution_search = neighbours(max_neighbour_index)
        else:
            break

    # Check which solution is better
    solution_values = [solution_value_ilp,
                       solution_value_greedy,
                       solution_value_search]
    solutions = [solution_ilp,
                 solution_greedy,
                 solution_search]
    best_idx = solution_values.index(max(solution_values))
    solution_value = solution_values[best_idx]
    solution = solutions[best_idx]

    # prepare the solution in the specified output format
    output_data = str(solution_value) + '\n'
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
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/ks_4_0)')

from collections import Counter
from random import choice, shuffle
from ortools.linear_solver import pywraplp

DEBUG = 0

def solve_it(input_data):
    # parse the input
    lines = input_data.split('\n')

    first_line = lines[0].split()
    node_count = int(first_line[0])
    edge_count = int(first_line[1])

    edges = []
    for i in range(1, edge_count + 1):
        line = lines[i]
        parts = line.split()
        edges.append((int(parts[0]), int(parts[1])))

    if DEBUG >= 1:
        print(f"Numero de vertices = {node_count}")
        print(f"Numero de arestas = {edge_count}")

    if DEBUG >= 2:
        print("Arestas:")
        for edge in edges:
            print(edge)
        print()

    return ColoringNaive(node_count, edge_count, edges)


def ColoringNaive(node_count, edge_count, edges):
    # Modify this code to run your optimization algorithm
    solutions = list()

    neighbours = dict()
    lookahead = dict()
    for node in range(node_count):
        neighbours[node] = {x for e in edges if node in e
                            for x in e if x != node}
        lookahead[node] = {x for n in neighbours[node]
                           for e in edges if n in e
                           for x in e if x != node}
    print('Lookahead pronto')

    # trivial solution: every node has its own color
    solution_trivial = range(0, node_count)
    solutions.append(solution_trivial)

    nodes_degrees = Counter()
    for s, t in edges:
        nodes_degrees[s] += 1
        nodes_degrees[t] += 1

    def run_randomised_alg(alg, k):
        best_sol = None
        best_sol_value = float('inf')
        iterations_without_change = 0

        while iterations_without_change <= k:
            cur_sol = alg()
            cur_sol_value = max(cur_sol) + 1
            if cur_sol_value < best_sol_value:
                best_sol = cur_sol
                best_sol_value = cur_sol_value
                iterations_without_change = 0
            else:
                iterations_without_change += 1
        return best_sol

    # Greedy, sorted, descending, randomised colour choice
    def greedy_descending_randomised():
        solution_greedy = [-1]*node_count
        for node, _ in nodes_degrees.most_common():
            neighbours_colours = {solution_greedy[n]
                                  for n in neighbours[node]
                                  if solution_greedy[n] >= 0}
            used_colours = {c for c in solution_greedy if c >= 0}
            available_colours = used_colours - neighbours_colours

            if available_colours:
                solution_greedy[node] = choice(list(available_colours))
            else:
                solution_greedy[node] = max(solution_greedy) + 1
        return solution_greedy

    def greedy_full_randomised():
        solution_greedy = [-1]*node_count
        shuffled_nodes = list(range(node_count))
        shuffle(shuffled_nodes)
        for node in shuffled_nodes:
            neighbours_colours = {solution_greedy[n]
                                  for n in neighbours[node]
                                  if solution_greedy[n] >= 0}
            used_colours = {c for c in solution_greedy if c >= 0}
            available_colours = used_colours - neighbours_colours

            if available_colours:
                solution_greedy[node] = choice(list(available_colours))
            else:
                solution_greedy[node] = max(solution_greedy) + 1
        return solution_greedy

    def greedy_adaptative():
        solution_greedy = [-1]*node_count

        # Set colour of the node with highest degree
        highest_node = nodes_degrees.most_common(1)[0][0]
        solution_greedy[highest_node] = 0

        colours_available = {n:set() for n in neighbours[highest_node]}

        while -1 in solution_greedy:
            if colours_available:
                # Assign colour to those nodes with less colours available
                cur_node = min(colours_available,
                               key=lambda x: len(colours_available.get(x)))
                if colours_available[cur_node]:
                    solution_greedy[cur_node] = choice(list(colours_available[cur_node]))
                else:
                    # New colour
                    new_colour = max(solution_greedy) + 1
                    solution_greedy[cur_node] = new_colour
                    for n in colours_available:
                        colours_available[n].add(new_colour)

                del colours_available[cur_node] # This node has already been processed
            else:
                # Nodes with no colour assigned can have any existing colours
                cur_node = choice([i for i, x in enumerate(solution_greedy) if x == -1])
                colour = choice(list({x for x in solution_greedy if x != -1}))
                solution_greedy[cur_node] = colour

            # Update dictionary of colours available
            for n in neighbours[cur_node]:
                if solution_greedy[n] == -1:
                    # Only process those nodes without colouring
                    if n in colours_available:
                        if solution_greedy[cur_node] in colours_available[n]:
                            colours_available[n].remove(solution_greedy[cur_node])
                    else:
                        colours_available[n] = {x for x in range(max(solution_greedy))
                                                if x != solution_greedy[cur_node]}
        return solution_greedy

    def ilp():
        max_degree = nodes_degrees.most_common(1)[0][1]

        solver = pywraplp.Solver('Colouring', pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
        colouring_var = [[solver.IntVar(0, 1, 'x[{}][{}]'.format(n, k))
                          for k in range(max_degree + 1)]
                         for n in range(node_count)]
        used_colour_var = [solver.IntVar(0, 1, 'y[{}]'.format(k))
                           for k in range(max_degree + 1)]

        for i, colours in enumerate(colouring_var):
            ct = solver.Constraint(1, 1, 'c_ct[{}]'.format(i))
            for var in colours:
                ct.SetCoefficient(var, 1)

        for s, t in edges:
            for k in range(max_degree + 1):
                ct = solver.Constraint(0, 1, 'e_ct[{}][{}][{}]'.format(s, t, k))
                ct.SetCoefficient(colouring_var[s][k], 1)
                ct.SetCoefficient(colouring_var[t][k], 1)

        for n in range(node_count):
            for k in range(max_degree + 1):
                ct = solver.Constraint(0, 1, 'u_ct[{}][{}]'.format(n, k))
                ct.SetCoefficient(colouring_var[n][k], -1)
                ct.SetCoefficient(used_colour_var[k], 1)

        objective = solver.Objective()
        objective.SetMinimization()
        for i, _ in enumerate(used_colour_var):
            objective.SetCoefficient(used_colour_var[i], 1)

        solver.Solve()
        solution = [[f.solution_value() for f in c].index(1) for c in colouring_var]

        colours_names = dict()
        cur_colour = 0
        for s in range(node_count):
            if solution[s] not in colours_names:
                colours_names[solution[s]] = cur_colour
                cur_colour += 1
            solution[s] = colours_names[solution[s]]
        return solution

    def greedy_lookahead():
        def lookahead_choice(node, colours, cur_sol):
            if len(colours) == 1:
                return colours[0]

            lookahead_colours = Counter()
            for n in lookahead[node]:
                if cur_sol[n] != -1:
                    lookahead_colours[cur_sol[n]] += 1
            # Get colour most present in the lookhead (neighbours of neighbours) to colour the current node
            return max(colours, key=lookahead_colours.get)
        solution_greedy = [-1]*node_count
        highest_node = nodes_degrees.most_common(1)[0][0]
        solution_greedy[highest_node] = 0

        colours_available = {n:set() for n in neighbours[highest_node]}

        while -1 in solution_greedy:
            if colours_available:
                cur_node = min(colours_available,
                               key=lambda x: len(colours_available.get(x)))
                if colours_available[cur_node]:
                    solution_greedy[cur_node] = lookahead_choice(cur_node,
                                                                 list(colours_available[cur_node]),
                                                                 solution_greedy)
                else:
                    new_colour = max(solution_greedy) + 1
                    solution_greedy[cur_node] = new_colour
                    for n in colours_available:
                        colours_available[n].add(new_colour)

                del colours_available[cur_node]
            else:
                cur_node = choice([i for i, x in enumerate(solution_greedy) if x == -1])
                colour = lookahead_choice(cur_node,
                                          list({x for x in solution_greedy if x != -1}),
                                          solution_greedy)
                solution_greedy[cur_node] = colour

            for n in neighbours[cur_node]:
                if solution_greedy[n] == -1:
                    if n in colours_available:
                        if solution_greedy[cur_node] in colours_available[n]:
                            colours_available[n].remove(solution_greedy[cur_node])
                    else:
                        colours_available[n] = {x for x in range(max(solution_greedy))
                                                if x != solution_greedy[cur_node]}
        return solution_greedy

    print('Greedy 1')
    solution_greedy1 = run_randomised_alg(greedy_descending_randomised, 5)
    solutions.append(solution_greedy1)
    print('Greedy 2')
    solution_greedy2 = run_randomised_alg(greedy_full_randomised, 5)
    solutions.append(solution_greedy2)
    print('Greedy Adaptative')
    solution_greedy3 = run_randomised_alg(greedy_adaptative, 5)
    solutions.append(solution_greedy3)
    print('Greedy Lookhead')
    solution_greedy_lokahead = greedy_lookahead()
    solutions.append(solution_greedy_lokahead)

    solutions_count = [max(s) + 1 for s in solutions]
    best_solution_idx = solutions_count.index(min(solutions_count))
    solution = solutions[best_solution_idx]
    n_colors = solutions_count[best_solution_idx]

    # prepare the solution in the specified output format
    output_data = str(n_colors) + '\n'
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
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/gc_4_1)')

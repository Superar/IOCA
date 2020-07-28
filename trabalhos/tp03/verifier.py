import sys
import fileinput
from os import listdir
from os.path import isfile, join

DEBUG = 2


def parse_input(instance_path):
    instance_file = open(instance_path, "r")
    input_data = instance_file.read()
    instance_file.close()

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

    instance_sol_file = open(instance_path + ".sol", "r")
    solution_data = instance_sol_file.read()
    instance_sol_file.close()

    # parse the solution
    lines = solution_data.split('\n')

    firstLine = lines[0].split()
    sol_value = int(firstLine[0])
    secondLine = lines[1].split()

    # list_answer = list of ints indicating that node i is labeled with color list_answer[i]
    list_answer = list(map(int, secondLine))

    return check_feasibility(node_count, edge_count, edges, sol_value, list_answer)


def check_feasibility(node_count, edge_count, edges, sol_value, list_answer):

    if (node_count != len(list_answer)):
        print("Solucao nao eh valida, pois numero de vertices na solucao esta incorreto")
        exit(0)

    for node_color in list_answer:
        if node_color < 0:
            print(
                "Solucao nao eh valida, pois existem vertices rotulados com numero negativo")
            exit(0)

    if DEBUG >= 1:
        print(
            f"Mapeamento de quais cores estão associadas a cada vertice: {list_answer}")
        print(f"Quantidade indicada de cores: {sol_value}")
        print(f"Quantidade real de cores: {(max(list_answer)+1)}")
        print()

    if sol_value != (max(list_answer)+1):
        print(
            "Erro! Quantidade de cores utilizadas nao corresponde ao indicado na solucao.")
        exit(0)

    for (u, v) in edges:
        if list_answer[u] == list_answer[v]:
            print(
                f"Solucao nao eh valida, pois vertices {u} e {v} sao vizinhos e tem a mesma cor")
            exit(0)

    # it is a feasible solution =D
    return sol_value


if __name__ == '__main__':
    if len(sys.argv) > 1:
        instance_path = sys.argv[1].strip()
        print(instance_path)
    else:
        print('This verifier requires an input file and a test type (0-2).  Please select one from the data directory (i.e. python verifier.py ./data/gc_50_3 2)')
        exit(0)
    # instance_path = input().rstrip()

    sol_value = parse_input(instance_path)

    test_type = int(sys.argv[2])
    # test_type = int(input())

    if test_type == 0:
        print("Solucao eh valida.")
        exit(0)

    instance_list = ['gc_50_3', 'gc_70_7', 'gc_100_5',
                     'gc_250_9', 'gc_500_1', 'gc_1000_5']
    good_values = [8, 20, 21, 95, 18, 124]
    great_values = [6, 17, 16, 78, 15, 100]

    instance_path_list = instance_path.split('\\')
    instance_name = instance_path_list[-1]
    i = instance_list.index(instance_name)

    if test_type == 1:
        if (sol_value <= good_values[i]):
            print("Parabens! Solucao eh boa.")
        else:
            print(
                f"Solucao nao eh boa, pois nao atingiu valor {good_values[i]}")

    if test_type == 2:
        if (sol_value <= great_values[i]):
            print("Parabens! Solucao parece ser otima.")
        else:
            print(
                f"Solucao nao eh otima, pois nao atingiu valor {great_values[i]}")

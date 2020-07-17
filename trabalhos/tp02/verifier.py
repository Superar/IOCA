import sys
import fileinput
from os import listdir
from os.path import isfile, join
from collections import namedtuple
import math

Point = namedtuple("Point", ['x', 'y'])
Facility = namedtuple(
    "Facility", ['index', 'setup_cost', 'capacity', 'location'])
Customer = namedtuple("Customer", ['index', 'demand', 'location'])

DEBUG = 1


def length(point1, point2):
    return math.sqrt((point1.x - point2.x)**2 + (point1.y - point2.y)**2)


def parse_input(instance_path):
    instance_file = open(instance_path, "r")
    input_data = instance_file.read()
    instance_file.close()

    # parse the input
    lines = input_data.split('\n')

    parts = lines[0].split()
    facility_count = int(parts[0])
    customer_count = int(parts[1])

    facilities = []
    for i in range(1, facility_count+1):
        parts = lines[i].split()
        facilities.append(Facility(
            i-1, float(parts[0]), int(parts[1]), Point(float(parts[2]), float(parts[3]))))

    customers = []
    for i in range(facility_count+1, facility_count+1+customer_count):
        parts = lines[i].split()
        customers.append(Customer(
            i-1-facility_count, int(parts[0]), Point(float(parts[1]), float(parts[2]))))

    if DEBUG >= 1:
        print(f"Numero de possiveis instalacoes = {facility_count}")
        print(f"Numero de clientes = {customer_count}")

    if DEBUG >= 2:
        print("Instalacoes na ordem que foram lidas")
        for facility in facilities:
            print(facility)
        print()
        print("Clientes na ordem que foram lidos")
        for customer in customers:
            print(customer)
        print()

    instance_sol_file = open(instance_path + ".sol", "r")
    solution_data = instance_sol_file.read()
    instance_sol_file.close()

    # parse the solution
    lines = solution_data.split('\n')

    firstLine = lines[0].split()
    sol_value = float(firstLine[0])
    secondLine = lines[1].split()

    # list_answer = list of ints indicating that customer i is connected with facility list_answer[i]
    list_answer = list(map(int, secondLine))

    return check_feasibility(facility_count, facilities, customer_count, customers, sol_value, list_answer)


def check_feasibility(facility_count, facilities, customer_count, customers, sol_value, list_answer):

    if (customer_count != len(list_answer)):
        print("Solucao nao eh valida, pois numero de clientes na solucao esta incorreto")
        exit(0)

    for facility in list_answer:
        if facility < 0 or facility >= facility_count:
            print("Solucao nao eh valida, pois ha clientes nao sendo atendidos")
            exit(0)

    used = [0]*len(facilities)
    capacity_remaining = [f.capacity for f in facilities]
    cost_of_solution = 0

    customer_index = 0
    for facility_index in list_answer:
        used[facility_index] = 1
        capacity_remaining[facility_index] -= customers[customer_index].demand
        cost_of_solution += length(customers[customer_index].location,
                                   facilities[facility_index].location)
        customer_index += 1

        if capacity_remaining[facility_index] < 0:
            print(
                f"Solucao nao eh valida, pois a capacidade da instalacao {facility_index} é menor que a demanda de seus clientes")
            exit(0)

    # calculate the cost of the solution
    cost_of_solution += sum([f.setup_cost*used[f.index] for f in facilities])
    cost_of_solution = float(format(cost_of_solution, ".2f"))

    if DEBUG >= 1:
        print(
            f"Mapeamento de quais instalacoes estao atendendo cada cliente: {list_answer}")
        print(f"Custo indicado: {sol_value}")
        print(f"Custo real: {cost_of_solution}")
        print()

    # evitando erro de arredondamento numérico
    if (sol_value > cost_of_solution + 1.0) or (sol_value < cost_of_solution - 1.0):
        print("Erro! Custo nao corresponde ao indicado na solucao.")
        exit(0)

    # it is a feasible solution =D
    return cost_of_solution


if __name__ == '__main__':
    if len(sys.argv) > 1:
        instance_path = sys.argv[1].strip()
        print(instance_path)
    else:
        print('This verifier requires an input file and a test type (0-2).  Please select one from the data directory (i.e. python verifier.py ./data/fl_25_2 2)')
        exit(0)
    # instance_path = input().rstrip()

    sol_value = parse_input(instance_path)

    test_type = int(sys.argv[2])
    # test_type = int(input())

    if test_type == 0:
        print("Solucao eh valida.")
        exit(0)

    instance_list = ['fl_25_2', 'fl_50_6', 'fl_100_7', 'fl_100_1', 'fl_200_7', 'fl_500_7',
                     'fl_1000_2', 'fl_2000_2']
    good_values = [4000000, 4500000, 2050, 26000000, 5000000, 30000000,
                   10000000, 10000000]
    great_values = [3269822, 3732794, 1966, 22724066, 4711295, 27006099,
                    8879294, 7453531]

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

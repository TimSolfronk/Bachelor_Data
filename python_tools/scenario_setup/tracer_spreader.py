import argparse

tolerance = 1e-5

'''
A script that takes a python array of tracers (saved in a .txt file as example) and spreads the body coordinate (x) from the original one by {o} with {d} step-size
can be useful to find, where exactly the fault is transformed to for a non-planar fault
'''

def get_spread_tracer(line):
    if line.find("[") == -1:
        return [line]
    parts = line.split(",")
    x_coord = float((parts[0]).split("[")[-1])
    rest_line = ""
    for i in range(1,len(parts)):
        rest_line += "," + parts[i]
    
    return_lines = []
    spread_x_coord = x_coord - args.offset
    while spread_x_coord <= x_coord + args.offset + tolerance:
        return_lines.append("[" + str(round(spread_x_coord,3)) + rest_line)
        spread_x_coord += args.delta   

    return_lines.append("\n")
    return return_lines


parser = argparse.ArgumentParser(description='A script that takes a python array of tracers (saved in a .txt file as example) and spreads the body coordinate (x) from the original one by {o} with {d} step-size')
parser.add_argument("--f",               dest="file_name",      type=str,    required=True,     help="Name of the file with the current tracer coordinates")
parser.add_argument("--d",               dest="delta",          type=float,  required=True,     help="Step size of the tracer spreading")
parser.add_argument("--o",               dest="offset",         type=float,  required=True,     help="How far to spread in each direction")
args = parser.parse_args()

lines = []
with open(args.file_name, "r") as file:
    lines = file.readlines()

new_tracers = []
for line in lines:
    new_tracers += get_spread_tracer(line)

with open(args.file_name + "_spread.txt", "w") as file:
    file.writelines(new_tracers)
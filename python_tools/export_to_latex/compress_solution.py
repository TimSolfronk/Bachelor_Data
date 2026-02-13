import argparse
import sys
from pathlib import Path
import csv

# this file expects a csv file that was generated with "export_csv_to_latex-py"
parser = argparse.ArgumentParser(description='A script that compresses a csv file, by only keeping every s (e.g. third) row.')
parser.add_argument("--f",               dest="name",      type=str,   required=True,  help="Name of the file to be checked, does not need to contain '.csv'" )
parser.add_argument("--s",            dest="skip",          type=int, default=1,   help="only take every second row" )
args = parser.parse_args()

with open(args.name, newline="", encoding="utf-8") as f:
    csvData= list(csv.reader(f, delimiter=','))

# keep first row with names
output = []
output.append(csvData[0])
csvData = csvData[1:]

# go through data and skip rows
for i in range(len(csvData)):
    if i % args.skip != 0:
        continue
    output.append(csvData[i])

#save result in a new file with "compr{s}" in name
name = args.name.split("tracer")[0] + "compr" + str(args.skip) + "_tracer" + args.name.split("tracer")[1]

with open(name, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f,delimiter=",")
    writer.writerows(output)
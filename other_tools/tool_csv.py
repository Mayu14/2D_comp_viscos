# coding: utf-8
from pathlib import Path
import csv

def main():
    fname = "log.txt"
    res = []
    with open(fname, "r", encoding = "utf-8") as f:
        for line in f:
            res.append(Path(line.split(" ")[-1]).name.replace("\n", ""))
    
    print(len(res))
    with open("res.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(res)


if __name__ == '__main__':
    main()
import sys
import csv
import math
from collections import defaultdict

files = sys.argv[1:]
if len(files) == 0:
    print("No CSV files provided", file=sys.stderr)
    sys.exit(1)

sums = defaultdict(lambda: defaultdict(float))
counts = defaultdict(lambda: defaultdict(int))
non_numeric = defaultdict(dict)
header = None

for path in files:
    f = open(path, newline='')
    reader = csv.DictReader(f)
    if header is None:
        header = reader.fieldnames
    for row in reader:
        t = row.get('Time')
        if t is None:
            continue
        non_numeric[t]['Time'] = t
        for col in row:
            if col == 'Time':
                continue
            val = row[col]
            if not val:
                continue
            try:
                x = float(val)
                if math.isnan(x) or math.isinf(x):
                    if col not in non_numeric[t]:
                        non_numeric[t][col] = val
                else:
                    sums[t][col] += x
                    counts[t][col] += 1
            except:
                if col not in non_numeric[t]:
                    non_numeric[t][col] = val
    f.close()

all_times = set(sums.keys()) | set(non_numeric.keys())
try:
    times = sorted(all_times, key=lambda x: float(x))
except:
    times = sorted(all_times)

with open('average_simulation_output.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(header)
    for t in times:
        out_row = []
        for col in header:
            if col == 'Time':
                out_row.append(t)
            elif counts[t].get(col, 0) > 0:
                out_row.append(str(sums[t][col] / counts[t][col]))
            else:
                out_row.append(non_numeric[t].get(col, ''))
        w.writerow(out_row)

print('average_simulation_output.csv')

import csv

times = [7, 13, 19, 25]
virion = [1, 183, 574, 983]
dip = [0, 0, 10, 35]
both = [0, 1, 68, 713]
sus = [3551, 3368, 2900, 1821]
total = 3552

rows = []
for i in range(len(times)):
    t = times[i]
    v = virion[i]
    d = dip[i]
    b = both[i]
    s = sus[i]
    
    v_rate = v * 100.0 / total
    d_rate = d * 100.0 / total
    b_rate = b * 100.0 / total
    s_rate = s * 100.0 / total
    
    v_str = "%.2f" % v_rate
    d_str = "%.2f" % d_rate
    b_str = "%.2f" % b_rate
    s_str = "%.2f" % s_rate
    
    rows.append([t, v, d, b, s, total, v_str, d_str, b_str, s_str])

f = open("13_infection_counts_by_time.csv", 'w', newline='')
w = csv.writer(f)
w.writerow(['time', 'virion_counts', 'dip_counts', 'both_infected_counts',
            'susceptible_counts', 'total_cells', 'virion_rate', 'dip_rate',
            'both_infected_rate', 'susceptible_rate'])

for r in rows:
    w.writerow(r)

f.close()
print("Generated: 13_infection_counts_by_time.csv")
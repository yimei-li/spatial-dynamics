#!/usr/bin/env python3

import csv

def generate_latest_csv():
    data = [
        {
            'time': 7,
            'virion_counts': 1,
            'dip_counts': 0, 
            'both_infected_counts': 0,
            'susceptible_counts': 3551,
            'total_cells': 3552
        },
        {
            'time': 13,
            'virion_counts': 183,
            'dip_counts': 0,
            'both_infected_counts': 1, 
            'susceptible_counts': 3368,
            'total_cells': 3552
        },
        {
            'time': 19,
            'virion_counts': 574,
            'dip_counts': 10,
            'both_infected_counts': 68,
            'susceptible_counts': 2900,
            'total_cells': 3552
        },
        {
            'time': 25,
            'virion_counts': 983,
            'dip_counts': 35,
            'both_infected_counts': 713,
            'susceptible_counts': 1821,
            'total_cells': 3552
        }
    ]
    
    # Calculate rates
    for row in data:
        total = row['total_cells']
        row['virion_rate'] = f"{row['virion_counts']/total*100:.2f}"
        row['dip_rate'] = f"{row['dip_counts']/total*100:.2f}"
        row['both_infected_rate'] = f"{row['both_infected_counts']/total*100:.2f}"
        row['susceptible_rate'] = f"{row['susceptible_counts']/total*100:.2f}"
    
    # Write to CSV
    output_file = "13_infection_counts_by_time.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'virion_counts', 'dip_counts', 'both_infected_counts',
                        'susceptible_counts', 'total_cells', 'virion_rate', 'dip_rate',
                        'both_infected_rate', 'susceptible_rate'])
        
        for row in data:
            writer.writerow([
                row['time'],
                row['virion_counts'],
                row['dip_counts'], 
                row['both_infected_counts'],
                row['susceptible_counts'],
                row['total_cells'],
                row['virion_rate'],
                row['dip_rate'],
                row['both_infected_rate'],
                row['susceptible_rate']
            ])
    
    print(f"Generated: {output_file}")

if __name__ == "__main__":
    generate_latest_csv()
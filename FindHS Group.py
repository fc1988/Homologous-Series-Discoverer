import csv

input_filename = 'target.csv'
output_filename = 'output.csv'

mz_threshold = 1
rt_threshold = 0.05

result_combinations = []

grouped_data = {}
with open(input_filename, 'r', newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader) 
    for row in reader:
        mz, rt, kmd, hs_number = map(float, row)
        if hs_number not in grouped_data:
            grouped_data[hs_number] = []
        grouped_data[hs_number].append((mz, rt, kmd))

for hs_number, group_data in grouped_data.items():
    sorted_group_data = sorted(group_data, key=lambda x: x[0])
    group_combinations = []
    for i in range(len(sorted_group_data)):
        current_combination = [sorted_group_data[i]]
        current_mz, current_rt, current_kmd = sorted_group_data[i]
        for j in range(i+1, len(sorted_group_data)):
            next_mz, next_rt, next_kmd = sorted_group_data[j]
            if next_rt > current_rt and next_mz > current_mz:
                if abs(next_mz - current_mz) >= mz_threshold and abs(next_rt - current_rt) >= rt_threshold:
                    current_combination.append((next_mz, next_rt, next_kmd))
                    current_mz, current_rt, current_kmd = next_mz, next_rt, next_kmd
        if len(current_combination) >= 3:
            subset_exists = False
            for existing_combination in group_combinations:
                if set(current_combination).issubset(set(existing_combination[1:])):
                    subset_exists = True
                    break
            if not subset_exists:
                group_combinations.append([hs_number] + current_combination)
    
    result_combinations.extend(group_combinations)

with open(output_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Group ID', 'HS Number', 'mz', 'RT', 'KMD'])
    for idx, combination in enumerate(result_combinations, start=1):
        group_id = f'Group{idx}'
        hs_number = combination[0]
        for mz, rt, kmd in combination[1:]:
            writer.writerow([group_id, hs_number, mz, rt, kmd]) 

import csv

# 读取CSV文件
input_filename = 'target.csv'
output_filename = 'output.csv'

# 定义阈值
mz_threshold = 1
rt_threshold = 0.05

# 初始化结果列表
result_combinations = []

# 读取CSV文件并按照HS Number进行分组
grouped_data = {}
with open(input_filename, 'r', newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # 跳过标题行
    for row in reader:
        mz, rt, kmd, hs_number = map(float, row)
        if hs_number not in grouped_data:
            grouped_data[hs_number] = []
        grouped_data[hs_number].append((mz, rt, kmd))

# 处理分组数据，找出排列组合
for hs_number, group_data in grouped_data.items():
    sorted_group_data = sorted(group_data, key=lambda x: x[0])
    group_combinations = []
    for i in range(len(sorted_group_data)):
        current_combination = [sorted_group_data[i]]
        current_mz, current_rt, current_kmd = sorted_group_data[i]
        for j in range(i+1, len(sorted_group_data)):
            next_mz, next_rt, next_kmd = sorted_group_data[j]
            if next_rt > current_rt and next_mz > current_mz:
                # 检查 mz 和 rt 的差异是否超过阈值
                if abs(next_mz - current_mz) >= mz_threshold and abs(next_rt - current_rt) >= rt_threshold:
                    current_combination.append((next_mz, next_rt, next_kmd))
                    current_mz, current_rt, current_kmd = next_mz, next_rt, next_kmd
        # 检查组合长度并添加到结果列表
        if len(current_combination) >= 3:
            # 检查是否有更大的子集已经存在
            subset_exists = False
            for existing_combination in group_combinations:
                if set(current_combination).issubset(set(existing_combination[1:])):
                    subset_exists = True
                    break
            if not subset_exists:
                group_combinations.append([hs_number] + current_combination)
    
    # 添加到总的结果列表
    result_combinations.extend(group_combinations)

# 写入CSV文件
with open(output_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Group ID', 'HS Number', 'mz', 'RT', 'KMD'])  # 添加KMD列
    for idx, combination in enumerate(result_combinations, start=1):
        group_id = f'Group{idx}'
        hs_number = combination[0]
        for mz, rt, kmd in combination[1:]:
            writer.writerow([group_id, hs_number, mz, rt, kmd])  # 添加KMD数据到每行

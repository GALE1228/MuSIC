import os

def read_tsv_file(file_path):
    data = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.strip() == '':
                    continue
                fields = line.strip().split('\t')
                data.append(fields)
    except Exception as e:
        return None
    return data

def merge_and_calculate(variants_file, wt_file, output_file):
    variants_data = read_tsv_file(variants_file)
    wt_data = read_tsv_file(wt_file)
    
    if variants_data is None or wt_data is None:
        return False
    
    if len(variants_data) != len(wt_data):
        return False
    
    merged_data = []
    for i in range(len(variants_data)):
        if variants_data[i][0] != wt_data[i][0] or variants_data[i][1] != wt_data[i][1]:
            print(f"variants: {variants_data[i][0]}, {variants_data[i][1]}")
            print(f"wt: {wt_data[i][0]}, {wt_data[i][1]}")
            return False
        
        try:
            variant_value = float(variants_data[i][2])
            wt_value = float(wt_data[i][2])
            difference = variant_value - wt_value
        except (IndexError, ValueError) as e:
            return False
        
        merged_row = variants_data[i][:2]
        merged_row.append(str(variant_value))
        merged_row.append(str(wt_value))
        merged_row.append(str(difference))
        
        merged_data.append(merged_row)
    
    column_names = ["ID", "Position", "Variants_Value", "WT_Value", "Difference"]
    merged_data.insert(0, column_names)
    
    for i in range(len(merged_data)):
        merged_data[i].pop(1)
    
    merged_data[0] = ["ID", "Variants_Value", "WT_Value", "Difference"]
    
    try:
        with open(output_file, 'w') as file:
            for row in merged_data:
                file.write('\t'.join(row) + '\n')
        return True
    except Exception as e:
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='merge two predicted file result')
    parser.add_argument('--variants', required=True)
    parser.add_argument('--wt', required=True)
    parser.add_argument('--output')
    
    args = parser.parse_args()
    
    variants_file = args.variants
    wt_file = args.wt
    
    if args.output:
        output_file = args.output
    else:
        output_dir = os.path.dirname(variants_file)
        base_name = os.path.basename(variants_file).replace("_variants.inference", "_diff.inference")
        output_file = os.path.join(output_dir, base_name)
    
    merge_and_calculate(variants_file, wt_file, output_file)
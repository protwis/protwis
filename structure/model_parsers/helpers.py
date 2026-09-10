import csv

def csv_to_dict(file_path, header_in="row"):    
    with open(file_path, 'r') as f:
        csv_reader = csv.reader(f, delimiter=',', quotechar='"')
        if header_in == "row":
            header = next(csv_reader)
            values = next(csv_reader)
            if len(header) != len(values):
                raise ValueError(f"Header and values length mismatch in CSV file {file_path}.")
            return dict(zip(header, values))            
        elif header_in == "column":
            return_dict = {}
            for line in csv_reader:
                key, value = line[0], line[1]
                return_dict[key] = value
            return return_dict
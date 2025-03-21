import csv

input_file = "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Scripts/scripts_cellranger_8/SampleSheet_modified.csv"
output_file = "/data/scratch/LAB/enric/Proyecto_pitagoras/Analisis_pitagoras/Scripts/scripts_cellranger_8/SampleSheet_final.csv"

with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)

    data_section = False  # To identify when we're in the `[Data]` section
    for row in reader:
        # Check if we're in the `[Data]` section
        if len(row) > 0 and row[0].strip() == "Sample_ID":
            data_section = True

        # Truncate indices only in the `[Data]` section
        if data_section and len(row) > 4:
            # Truncate Index (column 4) and Index2 (column 5)
            row[3] = row[3][:8]  # Truncate first index to 8 nt
            row[4] = row[4][:8]  # Truncate second index to 8 nt
        writer.writerow(row)
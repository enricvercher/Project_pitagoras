#!/usr/bin/python

# Ruta con las lecturas de clonotipos exportadas de MiXCR
reads_dir = "/data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_16/Spatial_TCR/Processed_data/mixcr_3_output_hudson/reads"

# Archivo FASTQ de lecturas R1
fastqs = ["/data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_16/Spatial_TCR/Raw_data/SRR15052441_1.fastq.gz"]

# Lista de posiciones de tejido de Space Ranger
spatial_file = "/data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_16/Spatial_transcriptomics/Processed_data/spatial/tissue_positions_list.csv"

# Archivo de salida
counts_file = "/data/scratch/LAB/enric/zz_Spatial_Seq_Hudosn/Patient_16/Spatial_TCR/Patient16_output_counts_mixcr3_modified.tsv"

#import libraries
import glob
from collections import defaultdict
import unicodedata
import gzip
import io

#positions dictionary; spatial barcode is key, excluding the "-1"; value will be dictionary with count for each clone; this dictionary will be written to output
positions = {} 

#get list of positions on spatial transcriptomics slide
with open(spatial_file, "r") as f:
	for line in f:
		positions[line.split(",")[0].strip("-1")] = defaultdict(int)

#get reads for each clone
files = glob.glob(reads_dir + "/*.fastq.gz")

supporting_reads = defaultdict(list) #key is a clone, value is a list of supporting reads
read_to_clone = {} #key is a read name, value is the supporting clone; this is the inverse of supporting_reads dictionary
readnames = {} #all readnames; Yes if duplicated, No if not

#open the FASTQ files, add read names to supporting_reads and mark as duplicate if 2nd appearance
for file_ in files: #iterate through reads files in mixcr reads output directory
	clone = file_.split("_cln")[1].split(".")[0] #get clone number from filename
	with io.TextIOWrapper(io.BufferedReader(gzip.open(file_))) as f: #open zipped file
		for line in f:
			if line[0] == "@":
				read = str(line.strip("@").strip("\n")) #read name
				supporting_reads[clone].append(read) #add read name to dictionary containing supporting reads for each clone
				read_to_clone[read] = clone #add supported clone to dictionary of read names
				if readnames.get(read):
					readnames[read] = "Yes"
					read_to_clone.pop(read) #remove duplicate reads that support >1 clone
				else:
					readnames[read] = "No"

read_barcodes = {} #key is read, value is spatial barcode
used_umis = {} #umis previously detected go here; can be different UMIs per spot

save = False
for file_ in fastqs: #iterate through original read 1 FASTQ files
	with io.TextIOWrapper(io.BufferedReader(gzip.open(file_))) as f:
		for line in f:
			if save == True:
				barcode = str(line.strip("\n")[0:16]) #extract spatial barcode sequence
				umi = line.strip("\n")[16:28] #extract UMI sequence
				if not used_umis.get((barcode,umi)): #check to see if this spatial barcode/UMI combination has been previously observed; this is a rudimentary method of UMI correction
					used_umis[(barcode,umi)] = True #add spatial barcode/UMI combination to used_umis dictionary
					read = readname.replace("1:N:0", "2:N:0") 
					read_barcodes[read] = barcode #add spatial barcode to read name
					if positions.get(barcode) != None: #if the barcode matches a spatial barcode present on the slide
						clone = read_to_clone[read] #get the clone the read supports
						positions[barcode][clone] += 1 #add it to positions library
				save = False
			if line[0] == "@": #if name line of FASTQ
				readname = str(line.strip("@").strip("\n")) #get read name (is in the format 1:N:0)
				if readnames.get(readname.replace("1:N:0", "2:N:0")): #look for mate, not original read
					if readnames[readname.replace("1:N:0", "2:N:0")] == "No": #if the read is not duplicated
						save = True #if conditions above met, will proceed to line 59 above with seqence data
										
#save output  
with open(counts_file, "w") as f:
    # Escribir el encabezado con los clones como clone1, clone2, etc., comenzando con un tabulador
    f.write("\t" + "\t".join([f"clone{x}" for x in supporting_reads.keys()]) + "\n")  # Cambiar los nombres de los clones con formato 'clone1'

    
    # Escribir las filas con las posiciones y los valores correspondientes para cada clon
    for position in positions.keys():
        f.write(position + "\t")  # Escribir la posición (barcode espacial)
        for clone in supporting_reads.keys():
            f.write(str(positions[position][clone]) + "\t")  # Número de UMIs para cada clon en la posición
        f.write("\n")

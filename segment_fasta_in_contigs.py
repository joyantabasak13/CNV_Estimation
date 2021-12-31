import csv
dir="/content/drive/My Drive/"
inp_file_name = "reference.fa"
segment_file_name = "seg.csv"

file = open(dir+inp_file_name,"r") #second level input file which was output in first step
file1 = open(dir+inp_file_name+"_segmented"+".fa","w") #output file: reference.fa


seq=file.readline()
seq = file.read()
with open(dir+segment_file_name) as seg_csv_file:
    seg_reader = csv.reader(seg_csv_file, delimiter='\t')
    line_count = 0
    for row in seg_reader:
        if line_count == 0:
            #print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            #print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
            file1.write(">"+str(line_count)+"\n")
            #print(line_count)
            seg_start=int(row[2])-1
            seg_end=int(row[3])
            count=0
            #print(seq[seg_start:])
            for ch in seq[seg_start:seg_end]:
                file1.write(ch)
                count+=1
                if count%50 ==0:
                    file1.write("\n")
                    count = 0
            file1.write("\n")

            
            line_count += 1
    #print(f'Processed {line_count} lines.')


file.close()
file1.close()
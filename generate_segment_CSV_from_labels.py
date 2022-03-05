import csv
# (643169,644858) = 1690   2
# "ID","chrom","loc.start","loc.end","num.mark","seg.mean"
# "Sample.1",22,1,149150,149150,0.9998

class Label:
    def __init__(self, start, end, copy):
        self.start = start
        self.end = end
        self.copy = copy

    def print_me(self):
        print("%d %d %d" % (self.start, self.end, self.copy))


# def correct_labels(labels):
#     corrected_labels = []
#     for x in labels:
#         start_extra = int(x.start) // 50
#         end_extra = int(x.end) // 50
#         corrected_labels.append(Label(x.start - start_extra, x.end - end_extra, x.copy))
#     return corrected_labels


f_in_label = '/home/joyanta/Documents/Projects/CNV/DataProcessing/data/label_with_del_400_1.txt'
f_out_segments = '/home/joyanta/Documents/Projects/CNV/DataProcessing/data/ForceBreakpoint_segmented_result_Win_50_RCutoff_1000_LR_Alpha=0.0005_.csv'
f_in_segments = '/home/joyanta/Documents/Projects/CNV/DataProcessing/data/segmented_result_Win_50_RCutoff_1000_LR_Alpha=0.0005.csv'

labels = []
with open(f_in_label, 'r') as fin:
    for line_no, line in enumerate(fin, 1):
        words = [x.strip() for x in line.split(',')]
        if len(words) > 0:
            start = words[0][1:]
            second_word = [x.strip() for x in words[1].split(' ')]
            end = second_word[0][:-1]
            copy = second_word[5]
            labels.append(Label(int(start), int(end), int(copy)))

# labels = correct_labels(labels)

segments = [1]
count = 0
for x in labels:
    if int(x.start) > segments[count]:
        segments.append(int(x.start))
        count = count + 1
    segments.append(int(x.end))
    count = count + 1
    print(f"start {x.start} and end {x.end} ")

existing_segments = []
with open(f_in_segments, 'r') as fin:
    for line_no, line in enumerate(fin, 1):
        words = [x.strip() for x in line.split(',')]
        if line_no > 1:
            existing_segments.append(int(words[3]))

col_names = ["ID","chrom","loc.start","loc.end","num.mark","seg.mean"]
with open(f_out_segments, 'w') as csvfile:
    # creating a csv writer object
    csv_writer = csv.writer(csvfile, delimiter=',')
    # writing the fields
    csv_contents = [col_names]
    for i in range(len(segments)-1):
        # writing the data rows --> "Sample.1",22,1,149150,149150,0.9998
        row = []
        row.append("Sample.1")
        row.append(22)
        row.append(segments[i])
        row.append(segments[i+1]-1)
        row.append(row[3] - row[2] + 1)
        row.append(1)
        print(row)
        csv_contents.append(row)
    row = []
    row.append("Sample.1")
    row.append(22)
    row.append(segments[-1])
    row.append(existing_segments[-1])
    row.append(row[3] - row[2] + 1)
    row.append(1)
    print(row)
    csv_contents.append(row)
    csv_writer.writerows(csv_contents)

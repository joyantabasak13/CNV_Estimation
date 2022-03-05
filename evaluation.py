import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


class Label:
    def __init__(self, start, end, copy):
        self.start = start
        self.end = end
        self.copy = copy

    def print_me(self):
        print("%d %d %d" % (self.start, self.end, self.copy))


f_in_copyno = '/home/joyanta/Documents/Projects/CNV/DataProcessing/data/copyNo_actual.txt'
f_in_label = '/home/joyanta/Documents/Projects/CNV/DataProcessing/data/label_with_del_400_1.txt'
f_in_segments = '/home/joyanta/Documents/Projects/CNV/DataProcessing/data/segmented_result_Win_50_RCutoff_1000_LR_Alpha=0.0005.csv'
contig_copyno = {}
labels = []
segments = []
with open(f_in_copyno, 'r') as fin:
    for line_no, line in enumerate(fin, 1):
        words = line.split()
        if len(words) > 0:
            contig_copyno[words[2]] = int(words[5])
    # print(contig_copyno)

with open(f_in_label, 'r') as fin:
    for line_no, line in enumerate(fin, 1):
        words = [x.strip() for x in line.split(',')]
        if len(words) > 0:
            start = words[0][1:]
            second_word = [x.strip() for x in words[1].split(' ')]
            end = second_word[0][:-1]
            copy = second_word[5]
            labels.append(Label(int(start), int(end), int(copy)))

with open(f_in_segments, 'r') as fin:
    for line_no, line in enumerate(fin, 1):
        words = [x.strip() for x in line.split(',')]
        if line_no > 1:
            segments.append(int(words[3]))

# print(segments)
# print(contig_copyno)
print(len(segments))
print(len(contig_copyno))


def correct_labels(labels):
    corrected_labels = []
    for x in labels:
        start_extra = int(x.start) // 50
        end_extra = int(x.end) // 50
        corrected_labels.append(Label(x.start - start_extra, x.end - end_extra, x.copy))
    return corrected_labels


def correct_segments(segments):
    corrected_segments = []
    for x in segments:
        extra = int(x) // 50
        corrected_segments.append(int(x) - extra)
    return corrected_segments


#corrected_labels = correct_labels(labels)
corrected_labels = labels
segments = correct_segments(segments)
# for i in range(len(labels)):
#     labels[i].print_me()
#     corrected_labels[i].print_me()

# Actual segments
x1 = list(range(0, int(segments[-1]), 50))
print(len(x1))
y1 = [1] * len(x1)
y2 = [1] * len(x1)
prev = 0
for i in range(len(contig_copyno)):
    down_sampled = segments[i]//50
    y1[prev:down_sampled] = [contig_copyno[str(i)]] * (down_sampled - prev)
    prev = down_sampled
print("y1 done")

for x in corrected_labels:
    y2[x.start//50: x.end//50] = [x.copy] * ((x.end - x.start)//50)
print("y2 done")
print(len(x1))
print(len(y1))
print(len(y2))
# Bug occurs due to downsampling but effect is very low. only upto <20 points
if len(x1)>len(y2):
    for i in range(len(x1)-len(y2)):
        y2.append(1)

# x_seq = np.array(x1)
# y_pred = np.array(y1)
# y_actual = np.array(y2)
# print("CONVERSION DONE")
fig = plt.figure()
plt.xlabel("Sequence Length")
plt.ylabel("Copy Number")
plt.plot(x1, y2, label="Actual (Simulated)", alpha=.5)
print("1 plotted")
plt.plot(x1, y1, label="Predicted", alpha=.5)
print("2 plotted")
plt.legend()
plt.savefig("Downsampled_plotG1S1.png")
plt.show()

f_in = '/home/joyanta/Documents/Projects/CNV/Debug_demo/reference.fa'
f_out = '/home/joyanta/Documents/Projects/CNV/Debug_demo/demo_ref.fa'

with open(f_in, 'r') as fin, open(f_out, 'w') as fout:
    for lineno, line in enumerate(fin, 1):
        if lineno < 1000:
            fout.write(line)

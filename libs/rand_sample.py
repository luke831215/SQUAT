import sys
import random
import os

output_file, input_file, readsize, num_sample = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
num_sample = readsize if num_sample > readsize else num_sample
ids_file = '{0}.ids'.format(os.path.splitext(output_file)[0])

#print("sampling {0} out of {1} records".format(num_sample, readsize))
id_list = random.sample(range(readsize), num_sample)
id_list.sort()

with open(input_file) as infile:
    with open(output_file, 'w') as w1:
        with open(ids_file, 'w') as w2:
            cnt = 0
            crt_idx = 0
            for line in infile:
                if cnt == id_list[crt_idx]:
                    w1.write('@{}\n'.format(crt_idx))
                    #w1.write(line)
                    w1.write(infile.readline())
                    w1.write(infile.readline())
                    w1.write(infile.readline())
                    
                    identifier = line.strip()
                    w2.write('@{0}\t{1}\n'.format(crt_idx, identifier))

                    crt_idx += 1
                    if crt_idx == num_sample:
                        break
                else:
                    next(infile)
                    next(infile)
                    next(infile)
                cnt += 1

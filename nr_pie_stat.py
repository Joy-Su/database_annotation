import os
import re
import subprocess

import pandas as pd

def nr_pie_stat(input, output, limit_max, scriptpath):
    H = {}
    with open(input,'r') as f, open(output,'w') as out:
        for line in f:
            if line.startswith('#') or not line:
                continue
            A = line.split('\t')
            name = A[-1].rstrip()
            match = re.match(r'.*\[([^\]]+)\]$', name)
            if not match:
                continue
            else:
                full_name = match.group(1)
                H[full_name] = H.get(full_name, 0) + 1
        name = [i[0] for i in sorted(H.items(), key=lambda d: (d[1], d[0]), reverse=True)]
        count = {}
        name_sort = []
        for i in range(len(name)):
            if i<int(limit_max):
                count[name[i]] = H[name[i]]
                name_sort.append(name[i])
            else:
                count['other species'] = count.get('other species',0)+H[name[i]]
        name_sort.append('other species')
        out.write('Species_Name\tGene_Number\n')
        for k in name_sort:
            out.write(k+'\t'+str(count[k])+'\n')
    input_path = os.path.dirname(os.path.abspath(input))
    cmd = 'cd %s && /home/sujieyi/anaconda3/envs/r4.0/bin/Rscript %s -i %s' %(input_path, scriptpath ,output)
    subprocess.call(cmd, shell=True)
    if os.path.exists(os.path.join(input_path,'Rplots.pdf')):
        os.remove(os.path.join(input_path,'Rplots.pdf'))
if __name__ == '__main__':
    nr_pie_stat()




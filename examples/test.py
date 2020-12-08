from os import listdir
from os.path import isfile, join
my_path = './ms_maf'
onlyfiles = [f for f in listdir(my_path) if isfile(join(my_path, f))]
onlyfiles.remove('.DS_Store')

import csv

with open('./ms_maf/maf.tsv', 'wt') as out_file:
    tsv_writer = csv.writer(out_file, delimiter='\t')
    tsv_writer.writerow(['MAF'])
    for i in range(46):
        tsv_writer.writerow(['./examples/ms_maf/'+onlyfiles[i]])
    
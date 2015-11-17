'''
Execute vec_merge.cc for each vec file in a given folder 
and output a single vec file
'''

import os
import shutil
import argparse

def Merge(file_1, file_2, vec_file):
    os.system('../../build/utils/vec_merge -a {0} -b {1} -c{2}'.format(file_1, file_2, vec_file))



def MergeVecFiles(directory, vec_file_name):
    dir_list = os.listdir(directory)
    init_file =  '{0}/{1}'.format(directory, dir_list[0])

    MB1 = '{0}/MB1.vec'.format(directory)
    shutil.copy(init_file, MB1)
    dir_list.pop(0)

    
    MB2 = open(directory + "MB2.vec", "w")
    MB2.close()
    MB2 = os.path.join(directory, "MB2.vec")

    print directory
    for i,file in enumerate(dir_list):
        if i % 2:
            print os.path.join(directory, file)
            print MB2
            print MB1
            Merge(MB1, os.path.join(directory, file), MB2)
        else:
            Merge(MB2, os.path.join(directory,file), MB1)

    # if len(dir_list) % 2:
    #     Merge(MB2, MB1, vec_file_name)

    # else:
    #     Merge (MB1, MB2, vec_file_name)


def main():

    parser = argparse.ArgumentParser(description='Process Input')
    parser.add_argument('-d', '--dir', type=str, required=True, help='File Directory for image files must be specified')
    parser.add_argument('-v', '--vec', type=str, required=True, help='.vec file name must be specified')
    args = parser.parse_args()

    MergeVecFiles(args.dir, args.vec)

if __name__ == '__main__':
    main()
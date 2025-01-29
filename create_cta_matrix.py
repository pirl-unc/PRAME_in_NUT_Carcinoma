import argparse

from pprint import pprint



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True)
    parser.add_argument('-o', required=True)
    parser.add_argument('-l', required=True)

    return parser.parse_args()

def create_matrix(args):

    pep_matrix = {}

    pep_matrix['797'] = list('0'*(int(args.l)-1))

    with open(args.i, 'r') as ifo:
        for line in ifo.readlines():
           line = line.split()
           cl = str(line[0])
           pos = int(line[1])
           pep_len = len(line[2]) 
           print("cell line: {}".format(cl))
           if cl not in pep_matrix.keys():
               pep_matrix[cl] = list('0'*(int(args.l)-1))
           print("protein binary length: {}".format(len(pep_matrix[cl])))
           print("start pos (0-index): {}".format(pos-1))
           print("peptide length: {}".format(pep_len))
           print("current binary before application: {}".format(pep_matrix[cl][pos-1:pos-1+pep_len]))
           print("Current protein binary: {}".format(pep_matrix[cl]))
           for app_pos in range(pos-1, min(pos-1+pep_len, len(pep_matrix[cl]))):
               print("app pos: {}".format(app_pos))
               print(len(pep_matrix[cl]))
               if pep_matrix[cl][app_pos] == '0':
                   pep_matrix[cl][app_pos] = '1'

    return pep_matrix 

def print_output_matrix(pep_mtx, args):

   with open(args.o, 'w') as ofo:
       ofo.write("line	pos	peptide_status\n")
       for k,v in pep_mtx.items():
           for pos, status in enumerate(v):
               ofo.write("{}	{}	{}\n".format(k, pos+1, status))

def main():

    args = get_args()
    pep_mtx = create_matrix(args)
    print_output_matrix(pep_mtx, args)




if __name__ == '__main__':
    main()

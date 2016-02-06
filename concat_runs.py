

import argparse
import glob

#Concatenate runs into a single file with one header

#def process():

#    return None

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-in')
    
    args = parser.parse_args()

    for i, in_file in enumerate(glob.glob(getattr(args, 'in'))):
        
        with open(in_file, 'rb') as fh:

            with open('all_runs.tsv', 'a+') as out:
                if i == 0:                
                    for line in fh:
                        out.write(line)
                else:
                    fh.readline()
                    for line in fh:
                        out.write(line)



            


            

            





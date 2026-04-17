# Script transforms the DynamicTypeArray into an DynamicDatatypeArray
# for DynamicIntArray it needs to be called with -d int

# Input:
# -i inputfile
# -d datatype (use Java datatype, e.g. int, double)

import sys, getopt

def main(argv):
    # if no inputfile is given by parameter, default one is used
    inputfile = 'DynamicTypeArray.java'
    datatype = ''
        
    try:
        opts, args = getopt.getopt(argv,"hi:d:",["ifile=","dtype="])
    except getopt.GetoptError:
        print 'DynamicArrayTransform.py -i <inputfile> -d <datatype>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'DynamicArrayTransform.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-d", "--dtype"):
            datatype = arg

    # Uppercase of datatype is used for name
    className = 'Dynamic' + datatype[:1].upper() + datatype[1:] + 'Array'
    outputfile = className + '.java'
    
    with open(inputfile, "rt") as fin:
        with open(outputfile, "wt") as fout:
            for line in fin:
                line = line.replace('DynamicTypeArray', className)
                fout.write(line.replace('Type', datatype))

if __name__ == "__main__":
   main(sys.argv[1:])
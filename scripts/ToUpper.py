import glob

for fn in glob.glob('dataOut_June16/OrothologyGroupSeqs/*.fasta'):
    ifile = open(fn)
    lines = ifile.readlines()
    ifile.close()

    with open(fn,'wb') as ofile:
        for l in lines:
            if l.startswith('>'):
                n = l[1:].rstrip().upper()
                ofile.write(">{0}\n".format(n))
            else:
                ofile.write(l)
import csv, os, subprocess, re, sys
ifile = sys.argv[1]
with open(ifile,'r') as ifil:
    fil = ifil.readlines()
nLines = len(fil)
outdir = sys.argv[2]
if not os.path.isdir(outdir):
    os.mkdir(outdir)
outfile = outdir+'parsed'+ifile.split('/')[-1]
ofil = open(outfile, 'w')
wf = csv.writer(ofil, delimiter='\t')
for linenum in range(0, nLines):
    print linenum
    line = fil[linenum]
    if(line[0:18] == 'Group of orthologs'):
        orthoid = line.split('#')[1].split('.')[0]
        #print orthoid
        bestscore = line.split(' ')[-2]
        #print bestscore
        linenum = linenum+1
        line = fil[linenum]
        newlinestrs = line.split(':')
        celescore = newlinestrs[1].split(' ')[0]
        specscore = newlinestrs[-1][:-1]
        linenum = linenum+1
        line = fil[linenum]
        while(line[0:2] != '__' and linenum < nLines):
            print 'h'+str(linenum)
            if(line[0:9] == 'Bootstrap'):
                try:
                    linenum=linenum+1
                    line=fil[linenum]
                    continue                
                except:
                    continue
            linews = line.split('\t')
            #print(line)
            #print(line[0:10])
            #print(linews)
            if(linews[0][0] != ' '):
                celegene=linews[0]
                celecov=linews[1][:-1]
            if(linews[3][0] != ' '):
                spegene=linews[3]
                specov=linews[4][:-2]
            linenum= linenum+1
            wf.writerow([orthoid,celegene,spegene,bestscore,celescore,specscore,celecov,specov])
            try:
                line=fil[linenum]
            except:
                pass
        linenum=linenum+1        
ofil.close()

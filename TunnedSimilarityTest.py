'''
    This script compares two .asc files and computes the simmilarity metrics I
    and D for these grids.

    The functions normalize, writenorfiles, getI and getD are part of the
    NichePy library available at https://github.com/bastodian/NichePy.

    This script is designed to loop through a folder (the variable Feature_Class
    below, specified as an argument to this script) containing subfolders each
    with Maxent predictions in the form of .asc for multiple species (in this 
    case aust_30 and tel_30 hardcoded below). These files are compared to
    corresponding Maxent .asc predictions in a folder called default and I and
    D metrics are generated for each comparison. 
    
    Using this setup the simmilarity of tuned Maxent models can be compared to
    default models using the I and D metrics. 
'''


def normalize(pairlist):
    '''
        This function creates a dictionary, normalizedictionary, of normalized grid values given a pair of .asc files.

        The structure of normalizedictionary is:
        {path/to/normalized/output/directory:[[ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value], [list of normalized grid values]]}

        input: pairlist (a list of file names)
    '''
    normalizedictionary = {}

    for file in pairlist:
        newpath = os.path.join(os.path.split(file)[0],'normalized_grids', os.path.split(file)[1])
        if ".asc" in file:
            with open(file, 'r') as thisfile:
                mylines = thisfile.readlines()
                linesX = [line.strip() for line in mylines if line.strip()]
                linespartitionA = [[line for line in linesX[0:6]],[x.split(' ') for x in linesX[6:len(linesX) - 1]]]
                linespartition = [linespartitionA[0],[]]
                nodata = linespartition[0][5].strip('NODATA_value').strip()
                for x in linespartitionA[1]:
                    for i in x:
                        if i != '': linespartition[1].append(i)
                mylen = len(linespartition[1])
                for k in range(mylen):
                    if linespartition[1][k] != nodata: linespartition[1][k] = float(linespartition[1][k])
                sumX = sum(linespartition[1][j] for j in range(mylen) if linespartition[1][j] != nodata)
                nrX = [linespartition[0],[]]
                for y in range(mylen):
                    if linespartition[1][y] == nodata: nrX[1].insert(y, linespartition[1][y])
                    else: nrX[1].insert(y, linespartition[1][y]/sumX)
                testX = sum(nrX[1][k] for k in range(mylen) if nrX[1][k] != nodata)
                normalizedictionary[newpath] = [nrX[0], [str(x) for x in nrX[1]]]
    return normalizedictionary

def writenorfiles(normalizedictionary):
    '''
        This is a function for writing normalized .asc files.

        input: normalizedictionary (the output of normalize)
    '''
    for item in normalizedictionary.keys():
        mypath = os.path.split(item)[0]
        if not os.path.exists(mypath): os.mkdir(mypath)

    for path, mylist in normalizedictionary.items():
        if os.path.exists(path): continue
        else:
            with open(path, 'w') as writeme:
                writeme.write('\n'.join(mylist[0]))
                writeme.write('\n')
                for i in mylist[0]:
                    if "ncols" in i:
                        me = i.strip("ncols")
                        ncols = int(me.strip())
                    elif "nrows" in i:
                        you = i.strip("nrows")
                        nrows = int(you.strip())
                        for i in range(1,nrows):
                            writeme.write('  '.join(mylist[1][(i-1)*ncols:(i*ncols)]))
                            writeme.write('\n')
                    else: continue

def getI(normalized, pairlist):
    '''
        This is a function for computing the I metric using the formula: 1 - 1/2(Sum((sqrt(Px,i)-sqrt(Py,i))^2)
        Where Px,i and Py,i are the values from the normalized grids from species X and Y at cell i.

        input: model (the niche modleing algorithm), normalized (the output of normalize), pairlist (the pair for which I will be calculated)
    '''
    mydictionary = {}

    for path, mylist in normalized.items():
        nodata = [x.strip('NODATA_value').strip() for x in normalized[path][0] if 'NODATA_value' in x]
        mypath = path.split('normalized_grids')
        mydictionary[''.join(list((mypath[0].strip('/'), mypath[1])))] = [nodata, [x for x in normalized[path][1]]]
    
    #for i,j in  zip(mydictionary[pairlist[0]][1],mydictionary[pairlist[1]][1]):
    #    print(i, j)

    mysum = sum((sqrt(float(x))-sqrt(float(y)))**2 for x,y in zip(mydictionary[pairlist[0]][1],mydictionary[pairlist[1]][1]) if x != mydictionary[pairlist[0]][0][0] and y != mydictionary[pairlist[1]][0][0])
    I = 1.0 - 0.5 * mysum

    print(I)
    return I

def getD(normalized, pairlist):
    '''
        This is a function for computing the D metric using the formula: 1 - 1/2(Sum(|Px,i-Py,i|)
        Where Px,i and Py,i are the values from the normalized grids from species X and Y at cell i.

        input: model (the niche modleing algorithm), normalized (the output of normalize), pairlist (the pair for which D will be calculated)
    '''
    mydictionary = {}

    for path, mylist in normalized.items():
        nodata = [x.strip('NODATA_value').strip() for x in normalized[path][0] if 'NODATA_value' in x]
        mypath = path.split('normalized_grids')
        mydictionary[''.join(list((mypath[0].strip('/'), mypath[1])))] = [nodata, [x for x in normalized[path][1]]]
   
    mysum = sum(fabs(float(x) - float(y)) for x,y in zip(mydictionary[pairlist[0]][1],mydictionary[pairlist[1]][1]) if x != mydictionary[pairlist[0]][0][0] and y != mydictionary[pairlist[1]][0][0])
    D = 1.0 - 0.5 * mysum
    
    print(D)
    return D

def myloop(FC, sp):
    RM = os.listdir(FC)
    myfile = FC+'_I_values.csv'
    default_list = os.listdir('default/'+sp)
    
    with open(myfile,'a') as outfile:
        outfile.write(','.join(['pair','I','D\n']))
        for mult in RM:
            sp_list = os.listdir(os.path.join(FC,mult,sp))
            for a in sp_list:
                if a == 'normalized_grids':
	            print('skipping '+a)
	        else:
	            for b in default_list: 
                        if b == 'normalized_grids':
		            print('skipping '+b)			
		        else:
		            mytest = os.path.join(FC,mult,sp,a)
                            mydefault = os.path.join('default',sp,b)
                            print('normalizing '+a+' in '+FC+'/'+mult)
		            normalize_dic = normalize([mytest, mydefault])
                            writenorfiles(normalize_dic)
                            print('computing I and D')
                            myI = getI(normalize_dic,[mytest, mydefault]) 
                            myD = getD(normalize_dic,[mytest, mydefault])
		            mylist = [FC+'/'+mult+': '+a+' vs L/1.0: '+b, str(myI), str(myD)+'\n']
                            outfile.write(','.join(mylist))
import os
import sys
from math import sqrt, fabs

Feature_Class = sys.argv[1]
aust = 'aust_30'
tel = 'tel_30'

myloop(Feature_Class, aust)
myloop(Feature_Class, tel)


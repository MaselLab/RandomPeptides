


import csv



def scan_file(data, filename): #Scans from the TSV file: takes in a list and string
    
                               #The list is the intended storage for the data set
                               #and the string is the file name
    
    with open(filename) as f:
        
        reader = csv.reader(f, delimiter="\t")
        
        for line in reader:
            data.append(line)
            #print(line)
        '''
        for i in range(100):
            print(data[i])
        '''
    #check for format    
    #for j in range(4):
    #    print(data[j])


def dicVals(dataPull, pullList, pepBook): #Creates a list of the needed information about
                                          #each peptides being pulled from the data set.
    
                                          #[PEP#, [Sums each day], ISD]
    
                                          #Takes in the list of all peptides, a list of the
                                          #requested peptides, and a dictionary containing
                                          #the sums of the observed values for for all four
                                          #days.
    output2 = open('idGrid.txt','w')
    output3 = open('ISDVals.txt', 'w')
    output4 = open('pepNames.txt','w')
    for pepn in range(1, len(dataPull)):
        
        pepName = dataPull[pepn][1]
        ISD = dataPull[pepn][2]
        FitnessBin = dataPull[pepn][5]
        
        pullList.append([pepName, pepBook[pepName], ISD, FitnessBin])
    for i in range(len(pullList)):
        
        output2.write(str(pullList[i][0]) + ", " + 
                     str(pullList[i][1]) + ", " +
                     str(pullList[i][2]) + ", " +
                     str(pullList[i][3]) + '\n')
        output4.write(str(pullList[i][0])+ '\n')
        output3.write(str(pullList[i][2]) + '\n')
    output2.close()
    output3.close()
    output4.close()
        


def daySums(dataExp, daySumsReg, pepBook): #Creates a dictionary of the peptides and their
                                           #sums for future use, as well as a short list of
                                           #the total sums for each day.

                                           #Takes in the experimental data set (list) as
                                           #well as a list and dictionary for the produced
                                           #values.
       
    for pep in range(1, len(dataExp)):
        
        pepName = dataExp[pep][0]
        pepBook[pepName] = [0] * 4
        
        for day in range(4):
            
            for rep in range(5):
                
                if (day < 3):
                    
                    a = int(dataExp[pep][2+rep*2+day*10])
                    
                    pepBook[pepName][day] += a
                    daySumsReg[day] += a
                    
                else:
                    
                    a = int(dataExp[pep][31 + rep])
                    
                    pepBook[pepName][day] += a
                    daySumsReg[day] += a

    print(daySumsReg)

def print2file(pullList, daySumsReg):
    #Creates output file
    output = open('PeptideDataOut.txt', 'w')
    '''
    output.write('PEPTIDE ID:' + ' ' + 'DaySums ' + 'ISD' + ' ' + 'BinFit' +'\n')
    output.write('\n')
    '''
    for pep in range(len(pullList)):
        
        
        output.write(str(pullList[pep][1][0]) + ", " + 
                     str(pullList[pep][1][1]) + ", " +
                     str(pullList[pep][1][2]) + ", " +
                     str(pullList[pep][1][3]) + '\n')
        '''
        for i in range(4):
            
            if (i == 1):
                
                a = len(str(pullList[pep][1]))
                
                while a < 24:
                    
                    output.write(' ')
                    a += 1
                    
                a = 0
                
            output.write(str(pullList[pep][i]))
            output.write('\t')
        output.write('\n')
        '''
    output.write(str(daySumsReg[0]) + ", " + 
                     str(daySumsReg[1]) + ", " +
                     str(daySumsReg[2]) + ", " +
                     str(daySumsReg[3]) + '\n')
        
      

    output.close()

def main():

    
    #Assignment/Creation of lists:
    dataExp = []
    dataPull = []
    pullList = []
    daySumsReg = [0] * 4
    pepSum = [0] * 4
    
    #Assignment/Creation of dictionaries:
    pepBook = {}

    #Variable used for cross-checking outputs
    k = 0

    
    #Creates data sets (lists):
    scan_file(dataExp, 'exp.summary.dat')
    
    scan_file(dataPull, 'peptide_data_isd_fitness.tsv')


    #Main information processing functions:
    daySums(dataExp, daySumsReg, pepBook)
    dicVals(dataPull, pullList, pepBook)


    #Cross-check of info for first peptide
    #print(pepBook['PEPNR00000000002'])
    #print(dataPull[1][2])
    #print(k)

    print2file(pullList, daySumsReg)
        
        
    #print(pullList)
    '''
    secOut = open('idGrid.txt')
    for pep in range(len(pullList)):
        secOut.write(str(pullList[pep][0]))
        secOut.write(' ')
        secOut.write(str(pullList[pep][1]))
        secOut.write('\n')
    '''
    '''
    expOut = open('expOut.txt','w')
    for i in range(len(dataExp)):
        expOut.write(str(dataExp[i]))
        expOut.write('\n')
    '''
    # checks for totals
    '''
    a = 0
    b = 0
    c = 0
    d = 0
    for i in range(len(pullList)):
        print(str(pullList[i][1][0]) + ", " + 
                     str(pullList[i][1][1]) + ", " +
                     str(pullList[i][1][2]) + ", " +
                     str(pullList[i][1][3]) + '\n')
        a += pullList[i][1][0]
        b += pullList[i][1][1]
        c += pullList[i][1][2]
        d += pullList[i][1][3]
        print(str(a) + ' ' + str(b) + ' ' + str(c) + ' ' + str(d))
    print(len(dataExp[1]))
    '''
#Calls main
main()

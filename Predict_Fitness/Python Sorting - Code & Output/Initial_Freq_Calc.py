import csv



def scan_file(data, filename): #Scans from the TSV file: takes in a list and string
    
                               #The list is the intended storage for the data set
                               #and the string is the file name
    
    with open(filename) as f:
        
        reader = csv.reader(f, delimiter="\t")
        
        for line in reader:
            data.append(line)
            #print(line)
    
def daySums(dataExp, daySumsReg, pepBook): #Creates a dictionary of the peptides and their
                                           #sums for future use, as well as a short list of
                                           #the total sums for each day.

                                           #Takes in the experimental data set (list) as
                                           #well as a list and dictionary for the produced
                                           #values.
       
    for pep in range(1, len(dataExp)):
        
        for i in range(1, len(dataExp[pep])):
            daySumsReg += int(dataExp[pep][i])
            
    print(str(daySumsReg))
    
def pNot(fitness, dataSums, p0):
    totalDays = dataSums[len(dataSums)-1]
    print(str(totalDays))
    for i in range(len(dataSums)-1):
        pepTotal = sum(dataSums[i])
        denom = 0
        for j in range(len(totalDays)):
            denom += totalDays[j] * (1 + fitness[i][0])**(j+1)
        p0.append(pepTotal/denom)

def scan_file2(data, filename):
    file = open(filename,'r')
    lines = file.readlines()
    for i in range(len(lines)):
        line = lines[i].split(',')
        for i in range(len(line)):
            line[i] = float(line[i])
        data.append(line)
    file.close()
def main():

    
    #Assignment/Creation of lists:
    dataExp = []
    dataPull = []
    pullList = []
    fitness = []
    dataSums = []
    datGood = []
    daySumsReg = 0
    pepSum = [0] * 4
    p0 = []
    #Assignment/Creation of dictionaries:
    pepBook = {}

    
    #Creates data sets (lists):
    scan_file(dataExp, 'exp.summary.dat')
    scan_file2(dataSums, 'PeptideDataOut.txt')
    scan_file2(fitness, 'FitnessMath.txt')

    pNot(fitness, dataSums, p0)
    '''
    for i in range(len(fitness)):
        print(fitness[i][0])
    '''
    outPNot = open('outPNot.txt','w')
    
    for i in range(len(p0)):
        outPNot.write(str(p0[i]))
        outPNot.write('\n')
    
    #Sanity Checks (group 1):                                                 
    #print(len(dataExp[0]))
    #print(len(dataExp[1]))
    #print(dataExp[0][35])

    #Main information processing functions:
    #daySums(dataExp, daySumsReg, pepBook)
    
main()

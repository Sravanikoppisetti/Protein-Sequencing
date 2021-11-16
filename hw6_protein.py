"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file=open(filename)
    data=file.read().splitlines()
    data="".join(data)
    #print(data)
    return data

'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    initial_lst=[]
    for i in range(startIndex,len(dna),3):
       initial_lst.append(dna[i:i+3])
       if dna[i:i+3]=='TAG' or dna[i:i+3]=='TAA' or dna[i:i+3]=='TGA':
           break
    replacedletters=[]
    for i in initial_lst:
        x=i.replace("T","U")
        replacedletters.append(x)
    return replacedletters


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    CodonDict={}
    f=open(filename)
    a=json.load(f)
    for i,j in a.items():
        for k in j:
            CodonDict[k.replace("T","U")]=i
            #print(CodonDict)
    return CodonDict

'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein=[]
    if codons[0]=='AUG':
        protein.append('Start')
    for i in range(1,len(codons)):
        if codons[i] in codonD.keys():
            protein.append(codonD[codons[i]])
    return protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna_data=readFile(dnaFilename)
    codon_data=makeCodonDictionary(codonFilename)
    count=0
    synthesizeProteins=[]
    length=0
    while length<len(dna_data):
        if dna_data[length:length+3]=="ATG":
            rnd_data=dnaToRna(dna_data,length)
            protein=generateProtein(rnd_data,codon_data)
            synthesizeProteins.append(protein)
            length=length+3*len(rnd_data)
        else:
            length=length+1
            count+=1
    return synthesizeProteins


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    commonProteins=[]
    for i in proteinList1:
        for j in proteinList2:
            if i==j and i not in commonProteins:
                commonProteins.append(i)
    return commonProteins


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    combineProteins=[]
    for i in proteinList:
        for j in i:
            if i not in combineProteins:
                combineProteins.append(j)
    return combineProteins


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aminoAcid={}
    for i in aaList:
        if i not in aminoAcid:
            aminoAcid[i]=1
        else:
            aminoAcid[i]+=1
    return aminoAcid


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    human_list=combineProteins(proteinList1)
    elephant_lst=combineProteins(proteinList2)
    human_dict=aminoAcidDictionary(human_list)
    elephant_dict=aminoAcidDictionary(elephant_lst)

    AminoAcid=[]
    freq1={}
    freq2={}
    freq_diff=[]

    for i in human_dict:
        freq1[i]=human_dict[i]/len(human_list)
        if i not in AminoAcid and i!='Start' and i!='Stop':
            AminoAcid.append(i)

    for j in elephant_lst:
        freq2[j]=elephant_dict[j]/len(elephant_lst)
        if j not in AminoAcid and j!='Start' and j!='Stop':
            AminoAcid.append(j) 

    for k in AminoAcid:
        frequency1=0
        frequency2=0
        if k in freq1:
            frequency1=freq1[k]
        if k in freq2:
            frequency2=freq2[k]
        difference=frequency2-frequency1
        if difference < -cutoff or difference > cutoff:
            freq_diff.append([k,frequency1,frequency2])
    return freq_diff


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()
    #test.testReadFile()
    #test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    #test.testSynthesizeProteins()
    

    ## Uncomment these for Week 2 ##
    """
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    """
    # test.testCommonProteins()
    #test.testCombineProteins()
    # test.testAminoAcidDictionary()
    test.testFindAminoAcidDifferences()

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """

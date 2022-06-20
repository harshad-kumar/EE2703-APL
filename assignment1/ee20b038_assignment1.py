def rev_sentence(sentence): 
    words = sentence.split()  
    reverse_sentence = " ".join(list(reversed(words)))  
    return reverse_sentence 
from sys import argv, exit
CIRCUIT = '.circuit'
END = '.end'
if len(argv) != 2:
    print("Entire the file in right format as given below\n") ;print('\nUsage: %s <inputfile>' % argv[0])
    exit()
try:
    with open(argv[1]) as file:
        lines = file.readlines()        #stores the list of all
        start = -1; end = -2
        for line in lines:              # extracting circuit definition start and end lines
            if CIRCUIT == line[:len(CIRCUIT)]:
                start = lines.index(line)
            elif END == line[:len(END)]:
                end = lines.index(line)
                break
        if start >= end or start<0:                # validating circuit block
            print('Invalid circuit definition, use correct format')
            exit(0)
        inlist = []
        for line in lines[start+1:end]:
            for x in range((len(line))):
                if line[x] == '#':
                    line = line[:x]
                    break        
            inlist.append(line)
        lines.clear()    
        lines = inlist         
        for line in lines:
            if line[0] == 'R' or 'L' or 'C' or 'V' or 'I':
                if len(line.split()) != 4:
                    print("incorrect number of parameters in line",(line))    
                    exit(0)
                if line.split()[1].isalnum() != True or line.split()[2].isalnum() != True:
                    print("Invalid node designation, use only alphanumeric")  
                    exit(0)
            elif line[0] == 'E' or 'G':
                if len(line.split()) != 6:
                    print("incorrect number of parameters in line",(line))  
                    exit(0)
                if line.split()[1].isalnum() != True or line.split()[2].isalnum() != True or line.split()[2].isalnum() != True or line.split()[4].isalnum() != True:
                    print("Invalid node designation, use only alphanumeric") 
                    exit(0)     
            elif line[0] == 'H' or 'F':
                if len(line.split()) != 5:
                    print("incorrect number of parameters in line",(line))  
                    exit(0)
                if line.split()[1].isalnum() != True or line.split()[2].isalnum() != True:
                    print("Invalid node designation, use only alphanumeric")  
                    exit(0)
                if line.split()[3][0] != "V":
                    print("Enter proper voltage symbol in",line)
                    exit(0)         
            elif line[0] != 'R' or 'L' or 'C' or 'V' or 'I' or 'E' or 'G' or 'H' or 'F':
                print("Proper circuit elements are not entered.")
                exit(0)              
            outlist = [] 
        for line in lines:     
            line = rev_sentence(line) 
            outlist.append(line)
            #lines.clear()
            lines = reversed(outlist) 
        for line in lines:
                print(line)     
except IOError:
    print('Invalid file')
    exit()

#Gudivada Harshad Kumar EE20B038
import numpy as np
from sys import argv, exit
def exp_to_num(exp):
    a = exp.split('e')
    num = (float(a[0]))*(float(10**float(a[1])))
    return num
#As mentioned in point 1 create class to store different details of components to be used later
class Component:
    def __init__(self,component,name,node1,node2,value):
        self.component = component
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.value = value
if len(argv) != 2: #valiadating the input
    print("Entire the file in right format as given below\n") ;print('\nUsage: %s <inputfile>' % argv[0])
    exit()
CIRCUIT = ".circuit"
END = ".end"
AC = ".ac"
comp_list = [];p=0;pi = 3.14159265 #creating a list (comp_list) to store different class objects
try:
    with open(argv[1]) as file:
        lines = file.readlines()
        start = -1 ; end = -2 ; ac = -3 ; ac_is = False
        for line in lines:              # extracting circuit definition start and end lines
            if CIRCUIT == line[:len(CIRCUIT)]:
                start = lines.index(line)
            elif END == line[:len(END)]:
                end = lines.index(line)
            elif AC == line[:len(AC)]:  #Checking if AC source exists
                ac = lines.index(line)
                ac_is = True            
        if start >= end or start<0:                # validating circuit block
            print('Invalid circuit definition, use correct format')
            exit(0)
        try:
            if ac_is == True:
                a,b,freq = lines[ac].split("#")[0].split()
                if freq.isdigit() == True:
                    freq = 2*pi*float(freq)
                else:
                    freq = 2*pi*exp_to_num(freq)
            for line in lines[start+1:end]:
                name,node1,node2,*value = line.split("#")[0].split()

                if name[0] == 'R':
                    object = Component("Resistor",name,node1,node2,value)
                elif name[0] == 'C':
                    object = Component("Capacitor",name,node1,node2,value)
                elif name[0] == 'L':
                    object = Component("Inductor",name,node1,node2,value)
                elif name[0] == 'V':
                    object = Component("Voltage_Source",name,node1,node2,value)
                    p = p+1
                elif name[0] == 'I':
                    object = Component("Current_Source",name,node1,node2,value)   

                if len(object.value) == 1:
                    if object.value[0].isdigit() == True:
                        object.value = float(object.value[0])
                    else:
                        object.value = float(exp_to_num(object.value[0]))
                else:                                                     #Analysing power supply elements inputs
                    if ac_is == True and (object.name == "Voltage_Source" or "Current_Source"):
                        object.value = (float(object.value[1])/2)*complex(np.cos(float(object.value[2])),np.sin(float(object.value[2])))
                    elif ac_is == False and (object.name == "Voltage_Source" or "Current_Source"):
                        if object.value[1].isdigit() == True:
                            object.value = float(object.value[1])
                        else:
                            object.value = float(exp_to_num(object.value[1]))    
                comp_list.append(object)
                                         
        except IndexError:
            print("Enter proper values in netlist")
            exit()
#Linking nodes using dictionary
    Node = {}
    for object in comp_list:
        if object.node1 not in Node:
            if object.node1 == "GND":
                Node["n0"] = 0
            else:
                node = "n" + object.node1
                Node[node] = int(object.node1)
        if object.node2 not in Node:
            if object.node2 == "GND":
                Node["n0"] = 0
            else:
                node = "n" + object.node2
                Node[node] = int(object.node2)
    n = len(Node); q = 0            
# Create an array matrix and solve Mx = b,
    M = np.zeros(((n+p-1),(n+p-1)),dtype="complex_"); b = np.zeros(((n+p-1),1),dtype="complex_")     
    for object in comp_list:                                            #filling the matrix according to rules of nodal analysis
        if object.component == "Resistor":
            if object.node1 == "GND":
                M[int(object.node2)-1][int(object.node2)-1] += 1/object.value
            elif object.node2 == "GND":
                M[int(object.node1)-1][int(object.node1)-1] += 1/object.value 
            else:
                M[int(object.node1)-1][int(object.node1)-1] += 1/object.value
                M[int(object.node2)-1][int(object.node2)-1] += 1/object.value
                M[int(object.node1)-1][int(object.node2)-1] += -1/object.value
                M[int(object.node2)-1][int(object.node1)-1] += -1/object.value
        elif object.component == "Capacitor":
            if ac_is == True:
                Xc = -1/(float(object.value)*freq)
                object.value = complex(0,Xc)
                if object.node1 == "GND":
                    M[int(object.node2)-1][int(object.node2)-1] += 1/object.value
                elif object.node2 == "GND":
                    M[int(object.node1)-1][int(object.node1)-1] += 1/object.value
                else:
                    M[int(object.node1)-1][int(object.node1)-1] += 1/object.value
                    M[int(object.node2)-1][int(object.node2)-1] += 1/object.value
                    M[int(object.node1)-1][int(object.node2)-1] += -1/object.value
                    M[int(object.node2)-1][int(object.node1)-1] += -1/object.value
            elif ac_is == False:
                print("Enter energy storage devices like capacitor and inductor with proper ac voltage")
                exit()
        elif object.component == "Inductor":
            if ac_is == True:
                Xl = float((object.value)*freq)  
                object.value = complex(0,Xl)
                if object.node1 == "GND":
                    M[int(object.node2)-1][int(object.node2)-1] += 1/object.value
                elif object.node2 == "GND":
                    M[int(object.node1)-1][int(object.node1)-1] += 1/object.value
                else:
                    M[int(object.node1)-1][int(object.node1)-1] += 1/object.value
                    M[int(object.node2)-1][int(object.node2)-1] += 1/object.value
                    M[int(object.node1)-1][int(object.node2)-1] += -1/object.value
                    M[int(object.node2)-1][int(object.node1)-1] += -1/object.value
            elif ac_is == False:
                print("Enter energy storage devices like capacitor and inductor with proper ac voltage")
                exit()
        elif object.component == "Current_Source":
            if object.node1 == 'GND':
                b[int(object.node2)-1][0] -= object.value
            elif object.node2 == 'GND':
                b[int(object.node1)-1][0] += object.value
            else:
                b[int(object.node1)-1][0] += object.value
                b[int(object.node2)-1][0] -= object.value
        elif object.component == "Voltage_Source":
            if object.node1 == "GND":
                M[int(object.node2) - 1][n-1+q] +=1
                M[n-1+q][int(object.node2)-1] +=1
                b[n-1+q] += object.value
                q+=1
            elif object.node2 == "GND":
                M[int(object.node1) - 1][n-1+q] +=1
                M[n-1+q][int(object.node1)-1] +=1
                b[n-1+q] += object.value
                q+=1
            else:
                M[int(object.node2)-1][n-1+q] += 1
                M[int(object.node1)-1][n-1+q] -= 1    
                M[n-1+q][int(object.node2)-1] +=1
                M[n-1+q][int(object.node1)-1] -=1
                b[n-1+q] += object.value
                q+=1
    if np.linalg.det(M) == 0:
        print("The given circuit is not solvable")
        exit(0)
    else:
        V = np.linalg.solve(M,b)
        print(V,"\n")
        for i in range(n-1):
            print("Voltage at node",i+1,"=",V[i],"\n")
        for j in range(p):
            print("Current I in the voltage source",j+1,"=",V[j+n-1],"\n")                
except IOError:
    print('Invalid file')
    exit() 
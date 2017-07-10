from cvxopt import matrix, spdiag, solvers, sparse
import numpy as np

from test4 import nodes, lines, spotLoads, Z, transformers, slackbus

N = len(nodes) # number of nodes
L = len(lines) # number of lines

lineVariables = []
for i in range(0,len(lines)):
    line = lines[i]
    if i == 0:
        nPhases = len(line[4])
    for ph in range(0,nPhases):
        if line[4][ph] == 0:
            continue
        lineVariables.append([i,ph,'p'])
        lineVariables.append([i,ph,'q'])

nodeVariables = []
slackVariables = []
i = 0
for node in nodes:
    for ph in range(0,nPhases):
        if nodes[node][ph] == 0:
            continue
        nodeVariables.append([i,node,ph])
        if node == slackbus[0]:
            slackVariables.append(i)
        i += 1
print nodeVariables
print slackVariables
R = []
X = []

for line in lines:
    length = float(line[2])*0.00018939 # ft -> miles
    config = line[3]

    r = []
    x = []

    if config in transformers:
        for ph in range(0,nPhases):
            if line[4][ph] == 0:
                r.append(0.0)
                x.append(0.0)
            else:
                r.append(transformers[config][0])
                x.append(transformers[config][1])

    else:
        for ph in range(0,nPhases):
            r.append(length*Z[config][ph][0].real)
            x.append(length*Z[config][ph][0].imag)

    R.append(r)
    X.append(x)
    
M = len(nodeVariables)+len(lineVariables)

A = matrix(0.0,(M,M))
b = matrix(0.0,(M,1))

for variable in lineVariables:
    line = lines[variable[0]]
    ph = variable[1]
    kind = variable[2]

    nodei = line[0]
    nodej = line[1]

    A[variable[0],variable[0]] = 1.0 # line flow in question

    # now look for things connected to j
    for variable2 in lineVariables:
        if lines[variable2[0]][1] != nodej:
            continue

        if variable[1] != variable2[1]:
            continue

        if variable[2] != variable2[2]:
            continue

        # so variable 2 leaves from j and is the same phase and type as variable
        A[variable2[0],variable2[0]] = -1.0

    # lastly look for real loads at j

    if kind == 'p':
        offset = 0
    elif kind == 'q':
        offset = 1
        
    try:
        load = float(spotLoads[nodej][ph*2+offset])
    except:
        load = 0.0

    b[variable[0]] = load

    r = R[variable[0]][ph]
    x = R[variable[0]][ph]

    for node in nodeVariables:
        if node[2] != ph:
            continue

        if node[1] == nodei:
            Vi_index = node[0]
        elif node[1] == nodej:
            Vj_index = node[0]

    print Vj_index

    if kind == 'p':           
        A[len(lineVariables)+Vj_index,len(lineVariables)+Vi_index] = 1.0 # Vi   
        A[len(lineVariables)+Vj_index,len(lineVariables)+Vj_index] = 1.0 # Vj
        
        A[len(lineVariables)+Vj_index,variable[0]] = r # Pij  

    else:
        A[len(lineVariables)+Vj_index,variable[0]] = x # Qij
     
# and the slack bus
for V_index in slackVariables:
    A[len(lineVariables)+V_index,len(lineVariables)+V_index] = 1.0
    b[V_index] = float(slackbus[1])
    print V_index

   
sol = np.linalg.solve(A,b)
print sol
'''
linesPower = []
nodeVoltages = []

c = 0
for i in range(0,L):
    line = []
    for ii in range(0,2*nPhases):
        line.append(sol[c][0])
        c += 1
    linesPower.append(line)
for i in range(L,L+N):
    node = []
    for ii in range(0,nPhases):
        node.append(sol[c][0])
        c += 1
    nodeVoltages.append(node)

linesPower = matrix(linesPower)
print linesPower.T
nodeVoltages = matrix(nodeVoltages)
print nodeVoltages.T
'''



from cvxopt import matrix, spdiag, solvers, sparse
import numpy as np

#nodes = [650,646,645,632,633,634,611,684,671,692,975,652,680]
nodes = {650:0,646:1,645:2,632:3,633:4,634:5,611:6,684:7,671:8,692:9,975:10,
         652:11,680:12}

# Node A, Node B, Length (ft), Config.
lines = [[632,645,500,'603'],
         [632,633,500,'602'],
         [633,634,0,'XFM-1'],
         [645,646,300,'603'],
         [650,632,2000,'601'],
         [684,652,800,'607'],
         [632,671,2000,'601'],
         [671,684,300,'604'],
         [671,680,1000,'601'],
         [671,692,0,'Switch'],
         [684,611,300,'605'],
         [692,675,500,'606']]

# Node, Phase 1 (kW), Phase 1 (kVar), Phase 2 (kW), Phase 2 (kVar),
# Phase 3 (kW), Phase 3 (kVar)
spotLoads = {634:[160,110,120,90,120,90],
             645:[0,0,170,125,0,0],
             646:[0,0,230,132,0,0],
             652:[128,86,0,0,0,0],
             671:[385,220,385,220,385,220],
             675:[485,190,68,60,290,212],
             692:[0,0,0,0,170,151],
             611:[0,0,0,0,170,80]}

# Node A, Node B, Phase 1 (kW), Phase 1 (kVar), Phase 2 (kW), Phase 2 (kVar),
# Phase 3 (kW), Phase 3 (kVar)
distLoads = [[632,671,17,10,66,38,117,68]]

# in Ohnms per mile
Z = {'601':[[complex(0.3465,1.0179),complex(0.1560,0.5017),
             complex(0.1580,0.4236)],
            [complex(0.3375,1.0478),complex(0.1535,0.3849)],
            [complex(0.3414,1.0348)]],
     '602':[[complex(0.7526,1.1814),complex(0.1580,0.4236),
             complex(0.1560,0.5017)],
            [complex(0.7475,1.1983),complex(0.1535,0.3849)],
            [complex(0.7436,1.2112)]],
     '603':[[complex(0.0000,0.0000),complex(0.0000,0.0000),
             complex(0.0000,0.0000)],
            [complex(1.3294,1.3471),complex(0.2066,0.4591)],
            [complex(1.3238,1.3569)]],
     '604':[[complex(1.3238,1.3569),complex(0.0000,0.0000),
             complex(0.2066,0.4591)],
            [complex(0.0000,0.0000),complex(0.0000,0.0000)],
            [complex(1.3294,1.3471)]],
     '605':[[complex(0.0000,0.0000),complex(0.0000,0.0000),
             complex(0.0000,0.0000)],
            [complex(0.0000,0.0000),complex(0.0000,0.0000)],
            [complex(1.3292,1.3475)]],
     '606':[[complex(0.7982,0.4463),complex(0.3192,0.0328),
             complex(0.2849,-0.0143)],
            [complex(0.7891,0.4041),complex(0.3192,0.0328)],
            [complex(0.7982,0.4463)]],
     '607':[[complex(1.3425,0.5124),complex(0.0000,0.0000),
             complex(0.0000,0.0000)],
            [complex(0.0000,0.0000),complex(0.0000,0.0000)],
            [complex(0.0000,0.0000)]]}
          
N = len(nodes) # number of nodes
L = len(lines) # number of lines

# Acutally first I want to construct a Y matrix
'''
Y = matrix(complex(0.0,0.0),(N,N))

for i in range(0,L):
    line = lines(i)
    nodeA = line[0]
    nodeB = line[1]
    length = line[2]
    config = line[3]

    zPerMile = Z[config]

    yab = complex(0.0,0.0)
'''

# changed my mind, just going straight in with the constraints\
# first considering all of the Pij terms

# I need to determine my decision variables - ie. which Pij exist
# dumb version - assume 3L Pijs and 3L Qijs, in order of definition

#X = [0.0]*(6*L)

# starting the matrix of constraints
A = matrix(0.0,(6*L,6*L))
b = matrix(0.0,(6*L,1))

for l in range(0,L):
    line = lines[l]
    nodei = line[0]
    nodej = line[1]

    try:
        loads = spotLoads[nodej]
    except:
        loads = [0,0,0,0,0,0]

    for ph in range(0,3):
        A[6*l+ph,6*l+ph] = 1.0 # p
        A[6*l+ph+3,6*l+ph+3] = 1.0 # q

        # look for other things connected to j
        for l2 in range(0,L):
            if lines[l2][0] != nodej:
                continue

            # so we now know that line ii connects j to something
            A[6*l+ph,6*l2+ph] = -1.0
            A[6*l+ph+3,6*l2+ph+3] = -1.0

        b[6*l+ph] = float(loads[2*ph])
        b[6*l+ph+3] = float(loads[2*ph+1])

X = np.linalg.solve(A,b)

linesPower = []

c = 0
for i in range(0,L):
    line = []
    for ii in range(0,6):
        line.append(X[c][0])
        c += 1
    linesPower.append(line)

linesPower = matrix(linesPower)
print linesPower.T
    
    
    




import numpy as np
import scipy.sparse
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import subprocess
from matplotlib import cm
import pandas as pd
import math
import scipy.sparse.linalg as linalg
import time
import warnings
warnings.filterwarnings("error")

def readNodeT(file):
    nd = np.array([])
    ndFile = open(file+'.node')
    lines = ndFile.readlines()
    l = np.fromstring(lines[0], sep = ' ')
    nNodes = int(l[0])
    l = np.fromstring(lines[1], sep = ' ')
    nd = np.append(nd,l)
    for l in lines[2:]:
        if (l[0] != "#"):
            l = np.fromstring(l, sep = ' ')
            nd = np.vstack((nd,l))
    return (nd,nNodes)

def readElemT(file):
    el = np.array([],dtype = int)
    elFile = open(file+'.ele')
    lines = elFile.readlines()
    l = np.fromstring(lines[0], sep = ' ')
    nNodesElem = 2
    l = np.fromstring(lines[1], dtype = int, sep = ' ')
    el = np.append(el,l[1:3]-1)
    el = np.vstack((el,l[2:4]-1))
    el = np.vstack((el,l[1:4:2]-1))
    for l in lines[2:]:
        if (l[0] != "#"):
            l = np.fromstring(l,dtype = int, sep = ' ')
            el = np.vstack((el,l[1:3]-1))
            el = np.vstack((el,l[2:4]-1))
            el = np.vstack((el,l[1:4:2]-1))
    return (el,el.shape[0],2)

def fixDofBycoords(coord,dof):
        a

def plotElem(elem,node, lbName = False,alp = 1):
    xlim,ylim, _, _, _, _ = setlim(node)
    for cont,e in enumerate(elem):
        x1,x2 = 0,0
        path = mpath.Path
        path_data = [(path.MOVETO, [node[e[0]][0], node[e[0]][1]])]
        for j in range(e.size):
            path_data.append((path.LINETO,[node[e[j]][0], node[e[j]][1]]))
            x1 += node[e[j]][0]
            x2 += node[e[j]][1]
        x1/=e.size
        x2/=e.size
        path_data.append((path.CLOSEPOLY, [node[e[0]][0], node[e[0]][1]]))
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(path,color='grey', ec = 'black', alpha = alp)
        ax.add_patch(patch)
        if lbName:
            plt.annotate(str(cont), xy=(x1,x2), ha="center",bbox=dict(facecolor='white', edgecolor='black', alpha= 0.5, boxstyle='round,pad=0.1'))
        # plt.annotate(str(elemAngle[cont-1]), xy=(x1,x2), ha="center",bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.1'))
    
def plotElemDisp(elem,nodeNew,totalDisp):
    xlim,ylim, _, _, _, _ = setlim(nodeNew)
    norm = plt.Normalize(totalDisp.min(), totalDisp.max())
    colors = cm.rainbow(norm(totalDisp))
    for cont,e in enumerate(elem):
        x1,x2 = 0,0
        path = mpath.Path
        path_data = [(path.MOVETO, [nodeNew[e[0]][0], nodeNew[e[0]][1]])]
        for j in range(e.size):
            path_data.append((path.LINETO,[nodeNew[e[j]][0], nodeNew[e[j]][1]]))
            x1 += nodeNew[e[j]][0]
            x2 += nodeNew[e[j]][1]
        x1/=e.size
        x2/=e.size
        path_data.append((path.CLOSEPOLY, [nodeNew[e[0]][0], nodeNew[e[0]][1]]))
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(path,color=0.5*(colors[e[0]]+colors[e[1]]), ec =0.5*(colors[e[0]]+colors[e[1]]), alpha = 1)
        ax.add_patch(patch)

    fig.colorbar(cm.ScalarMappable(norm=norm,cmap='rainbow'), ax = ax, shrink = 0.7, location='left')
        

def plotNode(node,disp = [],totalDisp = [], pos = False, plotId = False,plotCoord = False):
    xlim,ylim, _, _, _, _ = setlim(node)
    if (not pos):
        for c,i in enumerate(node):
            plt.scatter(i[0],i[1],color = 'grey')
            if plotId:
                plt.annotate(str(c),xy =(i-0.02*(xlim+ylim)),color = 'black')
            if plotCoord:
                plt.annotate(str(f'({i[0]:.2f},{i[1]:.2f})'),xy =(i),color = 'black', fontsize = 8)
    else:
        nodeNew = node.copy()
        
        for i in range(nNodes):
            for j in range(nDofNode):
                nodeNew[i,j] += disp[nDofNode*i+j]
        norm = plt.Normalize(totalDisp.min(), totalDisp.max())
        colors = cm.rainbow(norm(totalDisp))
        for c,i in enumerate(nodeNew):        
            sct = plt.scatter(i[0],i[1],color = 'black')
            plt.annotate(str(c),xy =(i+0.01*(xlim+ylim)),color = 'black')
    setlim(node)


############################
#CALCULATE ELEMENT ROTATIONS
def get_angle(x1,y1,x2,y2):
    return math.degrees(math.atan2(y2-y1, x2-x1))

def getElemAngle(elem,node,nElem):
    elemAngle = np.zeros((nElem),dtype=int)

    for e in range(nElem):
        x1 = node[elem[e][0]][0] 
        y1 = node[elem[e][0]][1]
        x2 = node[elem[e][1]][0] 
        y2 = node[elem[e][1]][1]

        elemAngle[e] = get_angle(x1,y1,x2,y2)
    return(elemAngle)

# nd = pd.read_csv('nodes.txt', sep='\s+', header = None)
# nNodes = node.shape[0]
# node = nd.to_numpy()
# print(node)
# elem = np.loadtxt('elem.txt', dtype = int, delimiter=' ')
# nElem = elem.shape[0]
#=
def getStiffnessBeam(E,I,A,L):
    L2 = L*L
    L3 = L*L*L
    return np.array(
    [[E*A/L  ,0,0,-E*A/L,0,0],
    [0, 12*E*I/L3 , 6*E*I/L2 ,0,-12*E*I/L3, 6*E*I/L2 ],
    [0, 6*E*I/L2 , 4*E*I/L  ,0, -6*E*I/L2, 2*E*I/L  ],
    [ -E*A/L  ,0,0,  E*A/L  ,0,0],
    [0,-12*E*I/L3,-6*E*I/L2 ,0, 12*E*I/L3,-6*E*I/L2 ],
    [0, 6*E*I/L2 , 2*E*I/L  ,0,-6*E*I/L2, 4*E*I/L  ], 
    ])

def getStiffnessTruss(E,A,L):
    return E*A/L * np.array(
    [[1,0,-1,0],
     [0,0,0,0],
     [-1,0,1,0],
     [0,0,0,0],
    ])

def rotateT(mat,theta):
    s = np.sin(theta* np.pi / 180)
    c = np.cos(theta* np.pi / 180)
    rot = np.array(
        [[c,s,0,0],
        [-s,c,0,0],
        [0,0,c,s],
        [0,0,-s,c]]
    )
    return rot.T @ mat @ rot

def rotateBeam(mat,theta):
    s = np.sin(theta* np.pi / 180)
    c = np.cos(theta* np.pi / 180)
    rot = np.array(
       [[ c,s,0,0,0,0],
        [-s,c,0,0,0,0],
        [0,0,1,0,0,0],
        [0,0,0, c,s,0],
        [0,0,0,-s,c,0],
        [0,0,0,0,0,1]]
    )
    return rot.T @ mat @ rot

def mean(x1,x2):
    return [(x1[0]+x2[0])*0.5,(x1[1]+x2[1])*0.5]

def calcNorm(x,y):
    return math.sqrt(x*x+y*y)


def findNear(node,nNeighs = 1):
    nBinds = 0
    nearList = [[]]
    for i in range(nNodes):
        nearList.append([])
        distls = np.zeros(nNodes)
        for j in range(nNodes):
            if(i==j):
                distls[j] = 0
            distls[j] = math.dist(node[i],node[j])
            
        min = distls[0]
        ix = 0
        if(i == 0):
            min = distls[1]
        for j in range(nNodes):
            if(distls[j]<min and distls[j]!=0):
                ix = j
                min = distls[j]
        # print(f"min distance from node {i} is {min} from node {ix}")
        for j in range(nNodes):
            if(i<j):
                if(distls[j]<(1+nNeighs)*min and distls[j]!=0):
                    nearList[i].append(j)
                    nBinds+=1
    return nearList, nBinds


def bindnear(nearList):
    elem2 = np.array([[]], dtype = int)
    item = [0,nearList[0][0]]
    elem2 = np.append(elem2,item)
    for j in nearList[0][1:]:
        elem2 = np.vstack((elem2,[0,j]))
    for i in range(1,nNodes):
        for j in nearList[i]:
            if(i<j):
                elem2 = np.vstack((elem2,[i,j]))
    return elem2

def getElemBin(bitv,nearList,intersecting = []):
    elem = np.array([[0,0]],dtype = int)
    cont = 0
    for i,nd in enumerate(nearList):
        for j in nd:
            if (bitv[cont]):
                elem = np.vstack((elem,[i,j]))
            cont+=1
    elem = elem[1:]
    for i in intersecting:
        if bitv[i[0]]:
            for gene in i[1:]:
                bitv[gene] = False
    return elem



def plotElemSS(elem,nd,s,E):
    xlim,ylim, _, _, _, _ = setlim(nd)
    norm = plt.Normalize(s.min(), s.max())
    colors = cm.rainbow(norm(s))
    for cont,e in enumerate(elem):
        x1,x2 = 0,0
        path = mpath.Path
        path_data = [(path.MOVETO, [nd[e[0]][0], nd[e[0]][1]])]
        for j in range(e.size):
            path_data.append((path.LINETO,[nd[e[j]][0], nd[e[j]][1]]))
            x1 += nd[e[j]][0]
            x2 += nd[e[j]][1]
        x1/=e.size
        x2/=e.size
        path_data.append((path.CLOSEPOLY, [nd[e[0]][0], nd[e[0]][1]]))
        codes, verts = zip(*path_data)
        path = mpath.Path(verts, codes)
        patch = mpatches.PathPatch(path,color=colors[cont], ec =colors[cont])
        ax.add_patch(patch)
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap='rainbow'), ax = ax, shrink = 0.7, location='left', label = 'strain')
    s = s * E
    norm = plt.Normalize(s.min(), s.max())
    fig.colorbar(cm.ScalarMappable(norm=norm,cmap='rainbow'), ax = ax, shrink = 0.7, location='left', label = 'stress')
    
def plotConstraint(center, direction,xlim,ylim, xmax, xmin, ymax, ymin):

    Path = mpath.Path
    if(direction == 0):
        a = -0.025*(xlim)
        b = 0.025*(ylim)
        if((xmax - center[0])<(center[0]-xmin)):
            a=-a
            b= -b
        path_data = [
            (Path.MOVETO, [center[0], center[1]]),
            (Path.LINETO, [center[0]+a, center[1]-b]),
            (Path.LINETO, [center[0]+a, center[1]+b]),
            (Path.CLOSEPOLY, [center[0], center[1]]),
            (Path.MOVETO, [center[0]+1.25*a, center[1]-b]),
            (Path.LINETO, [center[0]+1.25*a, center[1]+b])
            ]
    else:
        a = 0.025*(xlim)
        b = 0.025*(ylim)
        path_data = [
            (Path.MOVETO, [center[0], center[1]]),
            (Path.LINETO, [center[0]-a, center[1]-b]),
            (Path.LINETO, [center[0]+a, center[1]-b]),
            (Path.CLOSEPOLY, [center[0], center[1]]),
            (Path.MOVETO, [center[0]-a, center[1]-b*1.25]),
            (Path.LINETO, [center[0]+a, center[1]-b*1.25])
            ]
    codes, verts = zip(*path_data)
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path,color='grey', ec = 'black')
    ax.add_patch(patch)

def plotFixedEnd(node,theta,n):

    s = 0.05*(xlim) *np.sin(theta*np.pi/180)
    c = 0.05*(ylim) * np.cos(theta*np.pi/180)

    Path = mpath.Path
    path_data =  [(Path.MOVETO, [node[n][0], node[n][1]]),
            (Path.LINETO, [node[n][0]-s, node[n][1]+c]),
            (Path.LINETO, [node[n][0]+s, node[n][1]-c]),
            (Path.CLOSEPOLY, [node[n][0], node[n][1]])
            ]

    codes, verts = zip(*path_data)
    path = mpath.Path(verts, codes)
    patch = mpatches.PathPatch(path,color='grey', ec = 'black', linewidth = 2)
    ax.add_patch(patch)


def plotLoad(center,direction,value,xlim,ylim):
    a=0.0
    b=1
    if(value < 0):
        a=-1
        b=0
    if(direction == 0):
        ar = mpatches.FancyArrow(center[0],center[1],dx=(a+b)*xlim*0.075 ,dy=0,length_includes_head=True,width = ylim*0.0125,head_length = xlim*0.025, color = 'orange', ec='black')
        plt.annotate(str(f'{value:.2f}'), xy=(center[0]+xlim*0.1*(a+b),center[1]+ylim*0.05*(-a+b)), ha="center")
    else:
        ar = mpatches.FancyArrow(center[0],center[1]-(a)*ylim*0.075,dx=0 ,dy=(a+b)*ylim*0.075,length_includes_head=True,width = 0.0125*xlim,head_length = 0.025*ylim, color = 'orange', ec='black')
        plt.annotate(str(f'{value:.2f}'), xy=(center[0],center[1]-ylim*0.1*(a-b)), ha="center")

    ax.add_patch(ar)
        

def plotDirichlet(nodeDofList, node, elemAngle = []):
    xlim,ylim, xmax, xmin, ymax, ymin = setlim(node)
    for i in range(nodeDofList.shape[0]):
        if(all(k == True for k in nodeDofList[i]) and model == 1):
            plotFixedEnd(node,elemAngle[i],i)
        else:
            for j in range(nodeDofList[i].size):
                
                if (nodeDofList[i][j]):
                    center = node[i]
                    plotConstraint(center,j,xlim,ylim, xmax, xmin, ymax, ymin)
            

def plotNewmann(nodeLoadList, node, docList = []):
    xlim,ylim, _, _, _, _ = setlim(node)
    docL = np.flip(np.append(docList,-1))
    cont = 0
    for i in range(nodeLoadList.size):
        while (cont+i == docL[-1]):
            docL = np.delete(docL,-1)
            cont +=1
        if(nodeLoadList[i]!=0):
            center = node[(i+cont)//nDofNode]
            plotLoad(center,(i+cont)%nDofNode,nodeLoadList[i],xlim,ylim)


def getL(nNodesElem,nodeL,node):
    x,y = np.array([]),np.array([])
    for i in range(nNodesElem):
        ix = nodeL[i]
        x = np.append(x,node[ix][0])
        y = np.append(y,node[ix][1])
    l = math.dist([x[0],y[0]],[x[1],y[1]])
    return l

def getObj(elem,node,nNodesElem):
    
    sumL = sum(list(map(lambda x: getL(nNodesElem,x,node),elem)))
    return sumL
    



#Stiffness matrix assemble
def globalstiffAssemble(node, elem,elemAngle,nElem,nNodes,nNodesElem,nDofNode,model,E,A,I):
    K = scipy.sparse.lil_matrix((nNodes*nDofNode,nNodes*nDofNode), dtype =np.double)
    if (model==0):
        lList = np.zeros(nElem)
        for i in range(nElem):
            nodeL = elem[i]
            theta = elemAngle[i]
            l = getL(nNodesElem,nodeL,node)
            lList[i] = l
            kEl = rotateT(getStiffnessTruss(E,A,l),theta)
            for j in range(nNodesElem):
                actualNode = nodeL[j]
                for k in range(nDofNode):
                    ix = nDofNode*actualNode+k
                    for icol in range(nNodesElem):
                        for jcol in range(nDofNode):
                            col = nDofNode*icol+jcol
                            K[ix,nDofNode*nodeL[icol]+jcol] += kEl[nDofNode*j+k,col] 
        return K, lList
    elif (model==1):
        for i in range(nElem):
            nodeL = elem[i]
            theta = elemAngle[i]
            l = getL(nNodesElem,nodeL,node)
            kEl = rotateBeam(getStiffnessBeam(E,I,A,l),theta)
            for j in range(nNodesElem):
                actualNode = nodeL[j]
                for k in range(nDofNode):
                    ix = nDofNode*actualNode+k
                    for icol in range(nNodesElem):
                        for jcol in range(nDofNode):
                            col = nDofNode*icol+jcol
                            K[ix,nDofNode*nodeL[icol]+jcol] += kEl[nDofNode*j+k,col]
        return K
    

def getBC(fixedDof,fixedNode,loadDof):
    doc = 0
    nodeDofList = np.zeros((nNodes,nDofNode),dtype =bool)
    nodeLoadList = np.zeros((nNodes*nDofNode), dtype = np.double)
    
    for i in loadDof:
        nodeLoadList[nDofNode*i[0]+i[1]] = i[2]

    for i in fixedNode:
        for j in range(nodeDofList[i].size):
            if (not nodeDofList[i,j]):
                nodeDofList[i,j] = 1
                doc +=1
            
    for i in fixedDof:
        if (not nodeDofList[i[0],i[1]]):
            nodeDofList[i[0],i[1]] = 1
            doc +=1

    return nodeDofList,nodeLoadList,doc
    
def removeFreeNodes(node,elem,nodeLoadList):
    freenodeList = np.empty(0)
    for i in range(node.shape[0]):
        if i not in elem and nodeLoadList[nDofNode*i:nDofNode*i+1] == 0:
            freenodeList = np.append(freenodeList,[nDofNode*i,nDofNode*i+1])
    return freenodeList 

def assignBC(K,doc,nodeLoadList,nodeDofList, nNodes, nDofNode, freenodeList = []):
    cont = 0
    docList= np.zeros(doc,dtype = int)
    for i in freenodeList:
        if not nodeDofList[int(i//nDofNode),int(i%nDofNode)]:
            docList = np.append(docList,int(i))

    for i in range(nNodes):
        for j in range(nDofNode):
            if (nodeDofList[i,j]):
                dof = nDofNode*i+j-cont
                docList[cont] = nDofNode*i+j
                K= scipy.sparse.hstack([ K[:,0:dof] , K[:,dof+1:] ])
                K = K.tolil()
                K= scipy.sparse.vstack([ K[0:dof,:] , K[dof+1:,:] ]) 
                K = K.tolil()
                cont+=1
    
    docList =np.sort(docList)
    for i in docList[::-1]:
        nodeLoadList = np.delete(nodeLoadList,i)


    K= K[K.getnnz(1)>0][:,K.getnnz(0)>0]

    return K,docList,nodeLoadList

def solve(K,nodeLoadList):
    return linalg.spsolve(K.tocsr(),nodeLoadList)


def calcStrainT(node, nodeNew,elem):
    strain = np.zeros(elem.shape[0])
    c=0
    for e in elem:
        v1n = nodeNew[e[0]]
        v2n = nodeNew[e[1]]
        v1o = node[e[0]]
        v2o = node[e[1]]
        l = math.dist(v1o,v2o)
        lnew = math.dist(v1n,v2n)
        strain[c] = (lnew-l)/l
        c+=1
    return strain

def calcStress(strain,E):
    return strain*E

def getDisp(disp,docList,nDofNode,nNodes):
    dispEx = np.zeros(nDofNode*nNodes)
    cont = 0 
    docL = np.flip(np.append(docList,-1))

    for i in disp:
        while (cont == docL[-1]):
            docL = np.delete(docL,-1)
            cont+=1
        dispEx[cont] = i
        cont+=1
    return dispEx

def scale(dispEx,fs):
    dispExmin = min(dispEx)
    scale= (1/abs(dispExmin)/10 * fs) if (abs(dispExmin) != np.nan and dispExmin!= 0) else 1

    return  dispEx * scale

def nodeLoads(nodeLoadList,nodeDofList,elem):
    for i,nd in enumerate(nodeDofList):
        for j,dof in enumerate(nodeDofList[i]):
            c = i*nDofNode+j
            if nodeLoadList[c] == 0 and not dof:
                if(i in (elem)):
                    nodeLoadList[c]+=0.1
    return nodeLoadList


def getTotalDisp(dispEx,nNodes,nDofNode):

    totalDisp = np.zeros(nNodes)
    for x in range(nNodes):
        totalDisp[x] = calcNorm(dispEx[nDofNode*x],dispEx[nDofNode*x+1])
    return totalDisp


def updatePos(node,dispEx,nNodes,nDofNode,dim):
    nodeNew = node.copy()
    for i in range(nNodes):
        for j in range(dim):
            d = dispEx[nDofNode*i+j]
            nodeNew[i,j] += d
    return nodeNew

def scaleNew (node,dispExScaled,nNodes,nDofNode,dim):
    nodeNewScaled = node.copy()
    for i in range(nNodes):
        for j in range(dim):
            nodeNewScaled[i,j] += dispExScaled[nDofNode*i+j]
    return nodeNewScaled

def setlim(n):
    xmin = min(n[:,0])
    xmax = max(n[:,0])
    ymin = min(n[:,1])
    ymax = max(n[:,1])

    xlen = abs(xmax-xmin)
    ylen = abs(ymax-ymin)
    if(xlen>ylen):                
        plt.xlim([xmin-xlen*0.2,xmax+xlen*0.2])
        plt.ylim([ymin-xlen*0.2,ymax+xlen*0.2])
        xlim = -xmin-xlen*0.2+xmax+xlen*0.2
        ylim = -ymin-xlen*0.2+ymax+xlen*0.2
    else:                
        plt.xlim([xmin-ylen*0.2,xmax+ylen*0.2])
        plt.ylim([ymin-ylen*0.2,ymax+ylen*0.2])
        xlim = xmax+xlen*0.2-xmin-xlen*0.2
        if xlim <0.3:
            xlim = 0.3
        
        ylim = ymin-ylen*0.2+ymax+ylen*0.2
        if ylim < 0.3:
            ylim = 0.3
    return xlim,ylim, xmax, xmin, ymax, ymin 

def optruss(config, elem,plot = False, plotFinal = False):

    E,A,mass,dim,node,nodeDofList,nodeLoadList_a,doc = config.unpack()
    global xlim,ylim, xmax, xmin, ymax, ymin 
    global nNodes, nElem, nDofNode, model
    global fig, ax
    nElem = elem.shape[0]
    model = 0
    nNodesElem = 2
    elemAngle = getElemAngle(elem,node,nElem)
    for i in range(node.shape[0]):
        if(i not in elem):
            if(nodeLoadList_a[i*nDofNode] != 0 or nodeLoadList_a[i*nDofNode+1]!= 0):
                # if plotFinal:
                #     plt.annotate('erro 1', xy = (1,1))
                #     plt.show()
                return False , False, False, False

    #Boundary Conditions

    K, lList = globalstiffAssemble(node,elem,elemAngle,nElem,nNodes,nNodesElem,nDofNode,model,E,A,-1)

    for i,e in enumerate(elem):
        gk = lList[i]*A*mass
        nodeLoadList_a[nDofNode*e[0]+1] -= gk*0.5
        nodeLoadList_a[nDofNode*e[1]+1] -= gk*0.5

    freeNodeList = removeFreeNodes(node,elem,nodeLoadList_a)
    K,docList,nodeLoadList = assignBC(K,doc,nodeLoadList_a,nodeDofList, nNodes, nDofNode,freeNodeList)

    try:
        nll = nodeLoadList.copy()
        nll[::2]+=0.0001
        solve(K,nll)
        nll2 = nodeLoadList.copy()
        nll2[1::2]-=0.0001
        disp = solve(K,nll2)
        disp = solve(K,nodeLoadList)
        dispEx = getDisp(disp,docList,nDofNode,nNodes)
    except:
        # if plotFinal:
        #     plt.annotate('erro 2', xy = (1,1))
        #     plt.show()
        return False , False, False, False

    totalDisp = getTotalDisp(dispEx,nNodes,nDofNode)
    nodeNew = updatePos(node,dispEx,nNodes,nDofNode,dim)
    strain = calcStrainT(node, nodeNew,elem)
    stress = calcStress(strain,E)
    if(plot):
        fig, ax = plt.subplots()
        plotNewmann(nodeLoadList,node,docList)
        plotDirichlet(nodeDofList,node,elemAngle)
        plotElem(elem,node)
        plt.show()
        if(plotFinal):
            dispExScaled = scale(dispEx,5) 
            nodeNewScaled = scaleNew(node,dispExScaled,nNodes,nDofNode,dim)
            fig, ax = plt.subplots()
            plotElemSS(elem,nodeNewScaled,strain,E)
            setlim(nodeNewScaled)
            plt.title("Stresses")
            plt.show()
    return dispEx,totalDisp,stress,True

##############################
#FEM 
def fem(inputOpt,inputBC,outputOpt,plotOpt,E,A,I,sf,mod,dim,fileTriangle,node,elem,fixedDof,fixedNode,coordsFixedDof,coordsFixedNode,loadDof):
    global xlim,ylim, xmax, xmin, ymax, ymin 
    global fig, ax
    global nNodes, nElem, nDofNode, model
    
    model = mod
    
    fig, ax = plt.subplots()
    triangle = inputOpt['Input with Triangle lib']

    nNodes = node.shape[0]
    if(inputOpt['Element Binding by Proximity']):
        nearList,nBinds = findNear(node)
        elem = bindnear(nearList)
        print(nearList)
        print(nBinds)

    nElem = elem.shape[0]


    #######
    #Pre Process


    #DoF per Node
    dofNodeList = [
        2, #Truss: 2 DoF per node
        3 #Bernoulli Euler: 3 DoF per node
    ]

    nDofNode = dofNodeList[model]

    #Nodes per Element
    nodeElemList = [
        2, #Truss: 2 nodes per Element
        2 #Bernoulli Euler: 2 nodes per Element
    ]

    nNodesElem = nodeElemList[model]

    if(triangle):
        file = './input/'+fileTriangle
        nod, nNodes = readNodeT(file)
        node = nod[:,1:3]
        elem, nElem, nNodesElem = readElemT(file)

    
    xlim,ylim, xmax, xmin, ymax, ymin = setlim(node)

    nNodes = node.shape[0]
    

    #########
    #Boundary Conditions
        
    doc = 0
    nodeDofList = np.zeros((nNodes,nDofNode),dtype =bool)
    nodeLoadList = np.zeros((nNodes*nDofNode), dtype = np.double)

    for i in loadDof:
            nodeLoadList[nDofNode*i[0]+i[1]] = i[2]

    for i in fixedNode:
        for j in range(nodeDofList[i].size):
            if (not nodeDofList[i,j]):
                nodeDofList[i,j] = 1
                doc +=1
            
    for i in fixedDof:
        if (not nodeDofList[i[0],i[1]]):
            nodeDofList[i[0],i[1]] = 1
            doc +=1
        
    if (inputBC['Add beam joints']):
        if (model == 1):
            for i,nd in enumerate(nodeDofList):
                if(not nodeDofList[i,2]):
                    nodeDofList[i,2] = True
                    doc +=1

    elemAngle = getElemAngle(elem,node,nElem)

    if outputOpt['plot Nodes']:
        xlim,ylim, xmax, xmin, ymax, ymin  = setlim(node)
        plId,plCd = 0,0
        if(plotOpt['Plot Nodes with node id']): plId = 1
        if(plotOpt['Plot Nodes with node coords']): plCd = 1            
        plotNode(node,plotId = plId, plotCoord=plCd)
        plt.show() 
        fig, ax = plt.subplots()

    if outputOpt['plot Elements']:
        xlim,ylim, xmax, xmin, ymax, ymin  = setlim(node)
        if plotOpt['Plot Elements with Constraints']:
            plotDirichlet(nodeDofList,node,elemAngle)
        if plotOpt['Plot Elements with Loads']:
            plotNewmann(nodeLoadList,node)
        plotElem(elem,node,lbName = plotOpt['Plot Elements with element id'])
        plt.show()
        fig, ax = plt.subplots()
    if outputOpt['Solve']:
        
        K = globalstiffAssemble(node,elem,elemAngle,nElem,nNodes,nNodesElem,nDofNode,model,E,A,I)
        freeNodeList = removeFreeNodes(node,elem,nodeLoadList)
        K,docList,nodeLoadList = assignBC(K,doc,nodeLoadList,nodeDofList, nNodes, nDofNode,freeNodeList)
        disp = solve(K,nodeLoadList)
        dispEx = getDisp(disp,docList,nDofNode,nNodes)
        fs = 3
        
        dispExScaled = scale(dispEx,fs)
        totalDisp = getTotalDisp(dispEx,nNodes,nDofNode)        
        nodeNew = updatePos(node,dispEx,nNodes,nDofNode,dim)
        nodeNewScaled = scaleNew(node,dispExScaled,nNodes,nDofNode,dim)
        strain = calcStrainT(node, nodeNew,elem)
        stress = calcStress(strain,E)

    
        
        ############
        #strain/stress
        if outputOpt['plot Stress and Strain']:
            if plotOpt['Plot Post Processing with Loads']:
                plotNewmann(nodeLoadList,nodeNewScaled,docList)
            if plotOpt['Plot Post Processing with Constraints']:
                plotDirichlet(nodeDofList,node,elemAngle)
            if plotOpt['Plot Post Processing with Undeformed configuration']:
                plotElem(elem,node,False,0.1)
            plotElemSS(elem,nodeNewScaled,strain,E)
            setlim(nodeNewScaled)
            plt.title("Strain/Stress")
            plt.show()
            fig, ax = plt.subplots()
            

        ##############
        #Displacements
        if outputOpt['plot Total Displacement']:
            if plotOpt['Plot Post Processing with Loads']:
                plotNewmann(nodeLoadList,nodeNewScaled,docList)
            if plotOpt['Plot Post Processing with Constraints']:
                plotDirichlet(nodeDofList,node,elemAngle)
            if plotOpt['Plot Post Processing with Undeformed configuration']:
                plotElem(elem,node,False,0.1)
            plotElemDisp(elem,nodeNewScaled,totalDisp)
            setlim(nodeNewScaled)
            plt.title("Displacements")
            plt.show()

        if outputOpt['Output csv with TotalDisplacement']:
            df = pd.DataFrame([])
            for i in range(nDofNode):
                df[f'u[{i}]'] = (dispEx[i::nDofNode])
            df.to_csv("./output/TotalDisplacement.csv",sep=';')        
            maxdisplacement = abs(df).idxmax().to_numpy()
            print("Max Displacements:")
            for i in range(nDofNode):
                c = mean(node[elem[maxdisplacement[i],0]],node[elem[maxdisplacement[i],1]])
                print(f'{i}: el. {maxdisplacement[i]} : ({c[0]:.2f},{c[1]:.2f}) = {df.iloc[maxdisplacement[i],i]}')
                
        if outputOpt['Output csv with Stress']:
            df = pd.DataFrame([])
            if (model == 0):
                df[f's[{i}]'] = stress
            if(model == 1):
                for i in range(nDofNode):
                    df[f'u[{i}]'] = (stress[i::nDofNode])
            df.to_csv("./output/Stress.csv",sep=';')
            maxstress = df.idxmax().to_numpy()
            minstress = df.idxmin().to_numpy()
            print("Max/Min Stresses:")
            if (model == 0):
                c = mean(node[elem[maxstress[0],0]],node[elem[maxstress[0],1]])
                cmin = mean(node[elem[minstress[0],0]],node[elem[minstress[0],1]])
                print(f'{i} - Max: el. {maxstress[0]} : ({c[0]:.2f},{c[1]:.2f}) = {df.iloc[maxstress[0],0]}')
                print(f'{i} - Min: el. {minstress[0]} : ({cmin[0]:.2f},{cmin[1]:.2f}) = {df.iloc[minstress[0],0]}')
            elif(model == 1):
                for i in range(nDofNode):
                    c = mean(node[elem[maxstress[i],0]],node[elem[maxstress[i],1]])
                    cmin = mean(node[elem[minstress[i],0]],node[elem[minstress[i],1]])
                    print(f'{i} - Max: el. {maxstress[i]} : ({c[0]:.2f},{c[1]:.2f}) = {df.iloc[maxstress[i],i]}')
                    print(f'{i} - Min: el. {minstress[i]} : ({cmin[0]:.2f},{cmin[1]:.2f}) = {df.iloc[minstress[i],i]}')
                
        if outputOpt['Output csv with Strain']:
            df = pd.DataFrame([])
            if (model == 0):
                df[f's[{i}]'] = strain
            if(model == 1):
                for i in range(nDofNode):
                    df[f'u[{i}]'] = (strain[i::nDofNode])
            df.to_csv("./output/Strain.csv",sep=';')
            maxstrain = df.idxmax().to_numpy()
            minstrain = df.idxmin().to_numpy()
            print("Max/Min Strains:")
            if (model == 0):
                c = mean(node[elem[maxstrain[0],0]],node[elem[maxstrain[0],1]])
                cmin = mean(node[elem[minstrain[0],0]],node[elem[minstrain[0],1]])
                print(f'{i} - Max: el. {maxstrain[0]} : ({c[0]:.2f},{c[1]:.2f}) = {df.iloc[maxstrain[0],0]}')
                print(f'{i} - Min: el. {minstrain[0]} : ({cmin[0]:.2f},{cmin[1]:.2f}) = {df.iloc[minstrain[0],0]}')
            elif(model == 1):
                for i in range(nDofNode):
                    c = mean(node[elem[maxstrain[i],0]],node[elem[maxstrain[i],1]])
                    cmin = mean(node[elem[minstrain[i],0]],node[elem[minstrain[i],1]])
                    print(f'{i} - Max: el. {maxstrain[i]} : ({c[0]:.2f},{c[1]:.2f}) = {df.iloc[maxstrain[i],i]}')
                    print(f'{i} - Min: el. {minstrain[i]} : ({cmin[0]:.2f},{cmin[1]:.2f}) = {df.iloc[minstrain[i],i]}')
                


from Bio.PDB import *
import numpy as np
from numpy import linalg as LA

parser = PDBParser()

#model：把每个原子的name，氨基酸名字，三维坐标 导入到dictionary里#
def Dicforstructure (structure):
    #for model in structure:
    #   print(f"model {model}")
    dic = {}
    #Define the number of Chain#
    model2 = structure[0]
    for chain in model2:
        #print(f"chain {chain}, Chain ID: {chain.id}")
        chainid1 = chain.id
    #Extract res number#
        chain_B = model2[chainid1]
        k = 1
        for res in chain_B:
            #print(f"Residue name: {res.resname}, number: {res.id[1]}")
            for atom in res:
                #print(f"{atom.name}, position:", atom.coord)
                dic [k] = [atom.name,res.resname, res.id[1], atom.coord]
                k = k+1
    #print("There are",k,"atoms in this protein")
    return dic

#把结构归零，找出结构中心归零
def structureadjustment_r(structure):
    sum = np.array([0,0,0])
    structure1 = {}
    for atom_coord_num in structure:
        sum = sum + structure[atom_coord_num]
    new = sum/(len(structure))
    #print(new)
    for atom_coord_num in structure:
        structure1[atom_coord_num] = structure[atom_coord_num] - new
    return structure1

#paper里的RMSD算法
def score(dic1,dic2, eigenvalue):
    residue = 0
    for k in dic1:
        residue_withousqrt = np.sum(np.square(dic1[k]) + np.square(dic2[k]))
        residue = (residue_withousqrt + residue)
    RMSD = ((residue - 2*eigenvalue)/len(dic1))**0.5
    return RMSD

#这个是用来测试的，这个计算方法是没有经过旋转处理的，只有归零
def score_nonrotation (dic1,dic2):
    residue = 0
    for k in dic1:
        residue_withousqrt = np.sum(np.square(dic1[k] - dic2[k]))
        #print(dic1[k], dic2[k], residue_withousqrt)
        residue = (residue_withousqrt + residue)
    RMSD = (residue/len(dic1)) ** 0.5
    return RMSD

#Rij的计算，paper里面有说
def eigenvaluemax_R11(structure1,structure2, i, j):
    R=0
    for k in structure1:
        vector1 = structure1[k]
        vector2 = structure2[k]
        R = R + vector1[i-1]*vector2[j-1]
    return R

def findeigenvalue(seq1,seq2):
    R11 = eigenvaluemax_R11(seq1,seq2, 1, 1)
    R12 = eigenvaluemax_R11(seq1,seq2, 1, 2)
    R13 = eigenvaluemax_R11(seq1,seq2, 1, 3)

    R21 = eigenvaluemax_R11(seq1,seq2, 2, 1)
    R22 = eigenvaluemax_R11(seq1,seq2, 2, 2)
    R23 = eigenvaluemax_R11(seq1,seq2, 2, 3)

    R31 = eigenvaluemax_R11(seq1,seq2, 3, 1)
    R32 = eigenvaluemax_R11(seq1,seq2, 3, 2)
    R33 = eigenvaluemax_R11(seq1,seq2, 3, 3)

    eigenmatrix1 = np.array([[R11+R22+R33,R23-R32,R31-R13,R12-R21]
                         ,[R23-R32,R11-R22-R33,R12+R21,R13+R31],
                          [R31-R13,R12+R21,R22-R11-R33,R23+R32]
                          ,[R12-R21,R13+R31,R23+R32,R33-R11-R22]])
    #print(eigenmatrix1)

    w,v = np.linalg.eig(eigenmatrix1)
    eigenvalues1 = w
    #print(eigenvalues1)
    max1 = eigenvalues1[0]
    for ele in range(0, len(eigenvalues1)):
        if (eigenvalues1[ele] > max1):
            max1 = eigenvalues1[ele]
    return max1



def alignment (structureP,structureG, i,j):
    if len(structureG) <= len(structureP):
        dicpre_seq1 = structureG
        dicgold_seq1 = structureP
    elif len(structureG) >= len(structureP):
        dicpre_seq1 = structureP
        dicgold_seq1 = structureG
    finaldicG_seq1 = {}
    finaldicP_seq1 = {}
    ffinaldicG_seq1 = {}
    ffinaldicP_seq1 = {}
    x_x = dicpre_seq1[1]
    y_y = dicgold_seq1[1]

    x = x_x[2]-1
    y = y_y[2] -1
    for eleP in dicpre_seq1:
        first2eleP = dicpre_seq1[eleP]
        for eleG in dicgold_seq1:
            first2eleG = dicgold_seq1[eleG]
            if first2eleG[1] == first2eleP[1]:
                if first2eleG[0] == first2eleP[0]:
                    if first2eleG[2] - y == first2eleP[2] - x:
                        finaldicG_seq1[eleP] = [first2eleG[3], first2eleG[2] -y]
                        finaldicP_seq1[eleP] = [first2eleP[3], first2eleP[2] - x]
                        continue

    if i != "/" and j !="/":
        for ele in finaldicP_seq1:
            a = finaldicP_seq1[ele]
            b = finaldicG_seq1[ele]
            if a[1] >= i and a[1] <= j:
                ffinaldicP_seq1[ele] = a[0]
                ffinaldicG_seq1[ele] = b[0]
    if i =="/" and j =="/":
        for ele in finaldicP_seq1:
            a = finaldicP_seq1[ele]
            finaldicP_seq1[ele] = a[0]
            b = finaldicG_seq1[ele]
            finaldicG_seq1[ele] = b[0]
        ffinaldicG_seq1 = finaldicG_seq1
        ffinaldicP_seq1 = finaldicP_seq1


    finaldicP_seq1_afterR = structureadjustment_r(ffinaldicP_seq1)
    finaldicG_seq1_afterR = structureadjustment_r(ffinaldicG_seq1)



    eigenvalue = findeigenvalue(finaldicP_seq1_afterR, finaldicG_seq1_afterR)
    result = score(finaldicP_seq1_afterR, finaldicG_seq1_afterR, eigenvalue)
    return result


#Input
structure_gold_seq = parser.get_structure("S1G", r"C:\Users\Zhen\Downloads\6w6wB.pdb")
structure_predict_seq = parser.get_structure("S1P", r"C:\Users\Zhen\Downloads\seq2p.pdb")
dicgold_seq = Dicforstructure(structure_gold_seq)
dicpre_seq = Dicforstructure(structure_predict_seq)

# ”/","/" is default to compare whole sequence. Users can replace them with their partial sequence beginning number i and ending number j.
print("RMSD:",alignment(dicgold_seq,dicpre_seq, "/","/"))


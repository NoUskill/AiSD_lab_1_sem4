import string
import heapq
from itertools import cycle

from Tools.scripts.objgraph import ignore


def Huffman_alg(s: str):
    class Huffman_root:
        def __init__(self, sim: str, freq: int):
            self.sim = sim
            self.freq = freq
            self.leftchild = None
            self.rightchild = None

        def __lt__(self, other):
            return self.freq < other.freq
        def __str__(self):
            return self.sim + ' ' + str(self.freq)
    def assem_huf_table(tables_huf_cod: dict, root: Huffman_root, link: str):
        if root.sim == None:
            assem_huf_table(tables_huf_cod, root.leftchild, link=link + '0')
            assem_huf_table(tables_huf_cod, root.rightchild, link=link + '1')
        else:
            tables_huf_cod[root.sim] = link
        return tables_huf_cod
    tables_freq={}
    for i in s:
        if i in tables_freq:
            tables_freq[i]=tables_freq[i]+1
        else:
            tables_freq[i]=1
    queue=[Huffman_root(key, tables_freq[key]) for key in tables_freq.keys()]
    heapq.heapify(queue)
    while len(queue)>1:
        temp_left=heapq.heappop(queue)
        temp_right = heapq.heappop(queue)
        mid=Huffman_root(None,temp_left.freq+temp_right.freq)
        mid.leftchild=temp_left
        mid.rightchild=temp_right
        heapq.heappush(queue,mid)
    tables_huf_cod={}
    return assem_huf_table(tables_huf_cod,queue[0],"")
def BWT(s_orig: str):
    def radix_sort(arr):
        for i in range(len(arr[0])-1,-1,-1):
            arr.sort(key=lambda x: x[i])
        return arr
    s=s_orig
    cycle_shifts=[]
    for i in range(len(s)):
        s=s[1:]+s[0]
        cycle_shifts.append(s)

    cycle_shifts=radix_sort(cycle_shifts)
    index=cycle_shifts.index(s_orig)
    BWT_string=""
    print(cycle_shifts)
    for i in range(len(s)):
        BWT_string+=cycle_shifts[i][-1]

    return BWT_string, index
def IBWT(s_BWT: str, index_orig: int):
    list_s_BWT=list(s_BWT)
    for i in range(len(list_s_BWT)):
        list_s_BWT[i]=list_s_BWT[i]+chr(i)
    list_s_BWT_sort = sorted(list_s_BWT)
    dict_s_BWT={}
    for i in range(len(list_s_BWT_sort)):
        dict_s_BWT[list_s_BWT[i]]=list_s_BWT_sort[i]
    temp_index = list_s_BWT[index_orig]
    orig_s = ""
    for _ in range(len(s_BWT)-1):
        orig_s += dict_s_BWT[temp_index][0]
        temp_index=dict_s_BWT[temp_index]
    orig_s+=list_s_BWT[index_orig][0]
    print(orig_s)
def MTF(s_orig:str):
    alphabet=sorted(set(s_orig))
    list_MTF=[]
    for s in s_orig:
        i=0
        for t in alphabet:
            if s==t:
                for j in range(i,0,-1):
                    alphabet[j],alphabet[j-1]=alphabet[j-1], alphabet[j]
                break
            i += 1
        list_MTF.append(i)
    return (list_MTF,alphabet)
def IMTF(list_MTF: list,alphabet:list):
    list_orig = []
    for i in list_MTF:
        list_orig.append(alphabet[int(i)])
        for j in range(int(i), 0, -1):
            alphabet[j], alphabet[j - 1] = alphabet[j - 1], alphabet[j]
    return list_orig
def RLE(s_orig:str):
    t=0
    s_RLE=[]
    for i in range(len(s_orig)-1):
        t+=1
        if s_orig[i]!=s_orig[i+1]:
            if t<3:
                s_RLE.extend(s_orig[i]*t)
            else:
                s_RLE.append(str(t))
                s_RLE.append(s_orig[i])
            t=0


    if((t+1)<3):
        s_RLE.extend(s_orig[-1] * (t+1))
    else:
        s_RLE.append(str(t + 1))
        s_RLE.append(s_orig[-1])
    return s_RLE
def IRLE(s_RLE:list):
    i=0
    s_orig=""
    while i < len(s_RLE):
        if(s_RLE[i].isdigit()):
            s_orig+=int(s_RLE[i])*s_RLE[i+1]
            i+=2
        else:
            s_orig +=s_RLE[i]
            i+=1
    print(s_orig)
def LZ77(s_orig:str):
    buffer=7
    s_dict = ""

    left_index=0
    right_index=buffer

    while(left_index != len(s_orig)):
        s_temp=s_orig[left_index:right_index]
        max_len_for_chek=min(len(s_temp), len(s_dict))
        while max_len_for_chek > 0 :
            s_temp_for_chek= s_temp[0:max_len_for_chek]
            right_index_dict=len(s_dict)
            left_index_dict= right_index_dict- max_len_for_chek
            while ( left_index_dict >=0):
                s_temp_dict=s_dict[left_index_dict:right_index_dict]
                if(s_temp_dict==s_temp_for_chek):
                    print(left_index-left_index_dict,max_len_for_chek, s_temp_for_chek, s_orig[left_index+len(s_temp_dict)])
                    left_index+=len(s_temp_dict)
                    right_index+=len(s_temp_dict)
                    s_dict +=s_temp_dict
                    max_len_for_chek=0
                    break
                left_index_dict-=1
                right_index_dict-=1
            max_len_for_chek-=1
        if max_len_for_chek ==0:
            print(1,0, s_orig[left_index])
        s_dict+=s_orig[left_index]
        left_index+=1
        if right_index+1 > len(s_orig):
            right_index= len(s_orig)
        else:
            right_index+=1



def read_byte_str(name: str, len_sim:int ):
    sim_str=""
    with open(name, "rb") as f:
        while byte := f.read(1):
            sim_str+=(chr(int(f'{ord(byte):08b}',2)))
    return sim_str


def cod_in_byte_str(s:str):
    return s.encode("utf-8", errors=ignore)



s=read_byte_str("test.exe",1)
print(s)

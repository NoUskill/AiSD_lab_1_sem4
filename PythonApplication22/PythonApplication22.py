import heapq
import pickle
from itertools import compress, cycle
from Tools.scripts.objgraph import ignore
class compres_HA:
    def __init__(self):
        self.s_orig=None
        self.Huf_table=None
    def compress_str(self,s_input:str,s_file_name_output:str):
        self.s_orig=s_input
        self.Huf_table=Huffman_alg(self.s_orig)
        s_cod=""
        for s in self.s_orig:
            s_cod+=self.Huf_table[s]
        str_cod=""
        n=8
        list_s=[s_cod[i:i + n] for i in range(0, len(s_cod), n)]
        for i in range(0,len(list_s)-1):
             str_cod+=chr(int(list_s[i],2))       
        if (len(list_s[-1]) != 8):
             num_len=8-len(list_s[-1])
             list_s[-1]+="0"*(num_len)
             str_cod+=chr(int(list_s[-1],2))
             str_cod=chr(int(num_len))+str_cod
        else:
            str_cod+=chr(int(list_s[-1],2))
            str_cod=chr(int(0))+str_cod      
        print_str_in_file(str_cod,s_file_name_output)
        return (self.Huf_table)
    def compress(self,s_file_name:str,s_file_name_output:str):
        self.s_orig=read_byte_str(s_file_name,8)
        self.Huf_table=Huffman_alg(self.s_orig)
        print("in_compress")
        s_cod=""
        for s in self.s_orig:
            s_cod+=self.Huf_table[s]
        str_cod=""
        list_s=[s_cod[i:i + n] for i in range(0, len(s_cod), n)]
        for i in range(0,len(list_s)-1):
             str_cod+=chr(int(list_s[i],2))       
        if (len(list_s[-1]) != 8):
             num_len=8-len(list_s[-1])
             list_s[-1]+="0"*(num_len)
             str_cod+=chr(int(list_s[-1],2))
             str_cod=chr(int(num_len))+str_cod
        else:
            str_cod+=chr(int(list_s[-1],2))
            str_cod=chr(int(0))+str_cod      
        print_str_in_file(str_cod,s_file_name_output)
    def decompres(self,s_file_name:str,s_file_name_input):
        s_cod=read_byte_str_bin(s_file_name_input,8)
        temp=""
        s_orig=""
        reversed_dict = {v: k for k, v in self.Huf_table.items()}
        num_null=int((s_cod[0:8]),2)
        if num_null==0:
            s_cod=s_cod[8:]
        else:
            s_cod=s_cod[8:-num_null]
        for i in s_cod:
            temp+=i
            if reversed_dict.get(temp) != None:
                s_orig+=reversed_dict.get(temp)
                temp=""
        print_str_in_file(s_orig,s_file_name)
    def decompres_str(self,s_file_name_input:str,Huf_table:dict):
        s_cod=read_byte_str_bin(s_file_name_input,8)
        temp=""
        s_orig=[]
        self.Huf_table=Huf_table
        reversed_dict = {v: k for k, v in self.Huf_table.items()}
        num_null=int((s_cod[0:8]),2)
        if num_null==0:
            s_cod=s_cod[8:]
        else:
            s_cod=s_cod[8:-num_null]
        for i in s_cod:
            temp+=i
            if reversed_dict.get(temp) != None:
                s_orig.append((reversed_dict.get(temp)))
                temp=""
        return s_orig
class compres_RLE:
    def compres(file_input:str,file_output:str):
        s=decode_and_read_from_file(file_input)
        encode_and_write_to_file(RLE(s,"˼"),file_output)
    def decompres(file_input:str,file_output:str):
        s=decode_and_read_from_file(file_input)
        encode_and_write_to_file(IRLE(s,"˼"),file_output)
class compres_BWT_RLE:
    def __init__(self):
        self.s=None
        self.index=[]
    def compres(self,file_input:str,file_output:str):
         #dic={}
         #for i in range(0,10240):
         #    dic[chr(i)]=0
         self.s=decode_and_read_from_file(file_input)
         l_orig,self.index=BWT(self.s)
         print("BWT-comp")
         s_orig="".join(l_orig)
         #s_temp=RLE(s_orig,"˼")
         #for s in s_temp:
         #    dic[s]+=1
         #for i in range(0,1024):
         #    print(chr(i),dic[chr(i)])
         encode_and_write_to_file(RLE(s_orig,"˼"),file_output)
         print("RLE-comp")
    def decompres(self,file_input:str,file_output:str):
        s_temp= decode_and_read_from_file(file_input)
        s_temp=IRLE(s_temp,"˼")
        print("IRLE-comp")
        s_orig=IBWT(s_temp,self.index)
        print("IBWT-comp")
        encode_and_write_to_file(s_orig,file_output)
class compres_BWT_MTF_HA:
    def __init__(self):
        self.s=None
        self.index=None
        self.alphabet=None
        self.Huf_table=None
    def compres(self,file_input:str,file_output:str):
        self.s=read_byte_str(file_input,8)
        l_orig,self.index=BWT(self.s)
        print("BWT-comp")
        s_BWT="".join(l_orig)   
        s_orig, self.alphabet=MTF(s_BWT)
        print("MTF-comp")
        temp=compres_HA()
        self.Huf_table,self.s=temp.compress_str(s_orig,file_output)
        print("HA-comp")
    def decompres(self,file_input,file_output):
        temp=compres_HA()
        s_MTF=temp.decompres_str(file_input,self.Huf_table)   
        if self.s == "".join(s_MTF):
            print("T")
        else:
            print("F")
        print("IHA-comp")
        s_o=IMTF(s_MTF,self.alphabet)
        print("IMTF-comp")
        s_orig=IBWT(s_o,self.index)
        print("IBWT-comp")
        print_str_in_file(s_orig,file_output)
class compres_BWT_MTF_RLE_HA:
    def __init__(self):
        self.s=None
        self.index=None
        self.alphabet=None
        self.Huf_table=None
    def compres(self,file_input:str,file_output:str):
        dic={}
        for i in range(0,512):
            dic[chr(i)]=0
        self.s=read_byte_str(file_input,8)
        l_orig,self.index=BWT(self.s)  
        print("BWT-comp")
        s_BWT="".join(l_orig)
        s_MTF, self.alphabet=MTF(s_BWT)
        print("MTF-comp")
        self.s=s_MTF
        for s in s_MTF:
            dic[s]+=1
        s_RLE=RLE(s_MTF,"㆖")
        for s in s_MTF:
            dic[s]+=1
        print("RLE-comp")
        temp=compres_HA()
        self.Huf_table=temp.compress_str(s_RLE,file_output)
        print("HA-comp")
        for i in range(0,512):
            print(chr(i),dic[chr(i)])
    def decompres(self,file_input,file_output):
        temp=compres_HA()
        l_RLE=temp.decompres_str(file_input,self.Huf_table)
        s_RLE="".join(l_RLE)
        print("IHA-comp")
        s_MTF=IRLE(s_RLE,"㆖")
        if s_MTF==self.s:
            print("T")
        else:
            print("=====")
            for i in range(0,len(self.s)):
                if s_MTF[i] != self.s[i]:
                    print(i)
                    print(s_MTF[i:i+8])
                    print(self.s[i:i+8])
                    break
            print(len(s_MTF),len(self.s))
            print("F")
            print("=====")
        print("IRLE-comp")
        l_MTF=list(s_MTF)
        s_o=IMTF(l_MTF,self.alphabet)
        print("IMTF-comp")
        s_orig=IBWT(s_o,self.index)
        print("IBWT-comp")
        print_str_in_file(s_orig,file_output)
class compres_LZ77_HA:
     def __init__(self):
         self.sim="㆖"
         self.s=None
         self.Huf_table=None
     def compres(self,file_input:str,file_output:str):
         self.s=read_byte_str(file_input,8)
         s_LZ77=LZ77(self.s,100,self.sim)
         print("LZ77-comp")
         temp=compres_HA()
         self.Huf_table=temp.compress_str(s_LZ77,file_output)
         print("HA-comp")
     def decompres(self,file_input:str,file_output:str):
         temp=compres_HA()
         l_LZ77=temp.decompres_str(file_input,self.Huf_table)
         s_LZ77="".join(l_LZ77)
         print("IHA-comp")
         s_orig=ILZ77(s_LZ77,self.sim)
         print_str_in_file(s_orig,file_output)
class compres_LZ77:
     def __init__(self):
         self.sim="˼"
         self.s=None
     def compres(self,file_input:str,file_output:str):
         self.s=decode_and_read_from_file(file_input)
         s_LZ77=LZ77(self.s,100,self.sim)
         print("LZ77-comp")
         encode_and_write_to_file(s_LZ77,file_output)
     def decompres(self,file_input:str,file_output:str):
         s_LZ77=decode_and_read_from_file(file_input)
         print("IHA-comp")
         s_orig=ILZ77(s_LZ77,self.sim)
         encode_and_write_to_file(s_orig,file_output)
class compres_LZ78_HA:
     def __init__(self):
         self.sim="㆖"
         self.s=None
         self.Huf_table=None
     def compres(self,file_input:str,file_output:str):
         self.s=read_byte_str(file_input,8)
         s_LZ77=LZ77(self.s,100,self.sim)
         print("LZ77-comp")
         temp=compres_HA()
         self.Huf_table=temp.compress_str(s_LZ77,file_output)
         print("HA-comp")
     def decompres(self,file_input:str,file_output:str):
         temp=compres_HA()
         l_LZ77=temp.decompres_str(file_input,self.Huf_table)
         s_LZ77="".join(l_LZ77)
         print("IHA-comp")
         s_orig=ILZ77(s_LZ77,self.sim)
         print_str_in_file(s_orig,file_output)
class compres_LZ78:
     def __init__(self):
         self.sim="˼"
         self.s=None
     def compres(self,file_input:str,file_output:str):
         self.s=decode_and_read_from_file(file_input)
         s_LZ78=LZ78(self.s,self.sim)
         print("LZ78-comp")
         encode_and_write_to_file(s_LZ78,file_output)
     def decompres(self,file_input:str,file_output:str):
         s_LZ78=decode_and_read_from_file(file_input)
         print("IHA-comp")
         s_orig=ILZ78(s_LZ78,self.sim)
         encode_and_write_to_file(s_orig,file_output)


def Huffman_alg(s):
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
def BWT(s_orig):
    n=10000
    list_s=[s_orig[i:i + n] for i in range(0, len(s_orig), n)]
    list_BWT_string=[]
    list_BWT_index=[]
    for l in list_s:
        s=l
        cycle_shifts=[]
        for i in range(len(s)):
            s=s[1:]+s[0]
            cycle_shifts.append(s)
        cycle_shifts=sorted(cycle_shifts)
        index=cycle_shifts.index(l)
        BWT_string=""
        for i in range(len(s)):
            BWT_string+=cycle_shifts[i][-1]
        list_BWT_string.append(BWT_string)
        list_BWT_index.append(index)
    return list_BWT_string, list_BWT_index
def IBWT(s_BWT: str, index_orig: list):
    n=10000
    list_s=[s_BWT[i:i + n] for i in range(0, len(s_BWT), n)]
    j=0
    s_result=""
    for t in list_s:
        list_s_BWT=list(t)
        for i in range(len(list_s_BWT)):
            list_s_BWT[i]=list_s_BWT[i]+chr(i)
        list_s_BWT_sort = sorted(list_s_BWT)
        dict_s_BWT={}
        for i in range(len(list_s_BWT_sort)):
            dict_s_BWT[list_s_BWT[i]]=list_s_BWT_sort[i]
        temp_index = list_s_BWT[index_orig[j]]
        orig_s = ""
        for _ in range(len(t)-1):
            orig_s += dict_s_BWT[temp_index][0]
            temp_index=dict_s_BWT[temp_index]
        orig_s+=list_s_BWT[index_orig[j]][0]
        j+=1
        s_result+= orig_s
    return(s_result)
def MTF(s_orig):
    alphabet=sorted(set(s_orig))
    alphabet_copy=alphabet.copy()
    s_MTF=""
    for s in s_orig:
        i=0
        for t in alphabet:
            if s==t:
                for j in range(i,0,-1):                 
                    alphabet[j],alphabet[j-1]=alphabet[j-1], alphabet[j]
                break
            i += 1
        s_MTF+=chr(i)
    return (s_MTF,alphabet_copy)
def IMTF(list_MTF:list,alphabet:list):
    list_orig = []
    for i in list_MTF:
        temp_i=int(ord(i))
        list_orig.append(alphabet[temp_i])
        for j in range(temp_i, 0, -1):
            alphabet[j], alphabet[j - 1] = alphabet[j - 1], alphabet[j]
    return list_orig
def RLE(s_orig,sim):
    t=0
    s_RLE=""
    for i in range(len(s_orig)-1):
        t+=1
        if s_orig[i]!=s_orig[i+1]:
            if t<5:
                s_RLE+=(s_orig[i]*t)
            else:
                s_RLE+=(sim+chr(t)+s_orig[i]+sim)
            t=0
    if((t+1)<5):
        s_RLE+=(s_orig[-1] * (t+1))
    else:
       s_RLE+=(sim+str(t + 1)+s_orig[-1]+sim)
    return (s_RLE)
def IRLE(s_RLE,sim):
    s_orig=""
    print("in IRLE")
    i=0
    while i<len(s_RLE):
        if s_RLE[i]==sim:
            temp_s=s_RLE[i:i+4]
            if temp_s[0]==temp_s[-1]==sim:
                s_orig+=temp_s[2]*ord(temp_s[1])
                i+=4
            else:
                s_orig+=s_RLE[i]
                i+=1
        else:
            s_orig+=s_RLE[i]
            i+=1
    return s_orig
def LZ77(s_orig,n,sim):
    def longest_prefix_from(Left, Right):
        LongestPrefixLength = 0
        LongestPrefixPos = -1
        while 1:
            PrefixLength = LongestPrefixLength+1
            if PrefixLength >= len(Right):
                break
            Prefix = Right[0: PrefixLength]
            PrefixPos = Left.find(Prefix)
            if PrefixPos == -1:
                break
            LongestPrefixLength = PrefixLength
            LongestPrefixPos = PrefixPos
        return (LongestPrefixLength, LongestPrefixPos)
    def codeBufferLZ77(Buffer):
        Result = ""
        CodePos = 0
        while CodePos < len(Buffer):
            Left = Buffer[0:CodePos]
            Right = Buffer[CodePos:]
            (PrefixLength, PrefixPos) = longest_prefix_from(Left, Right)
            if (PrefixLength == 0):
                Result+=(Buffer[CodePos])
                CodePos = CodePos+1
            else:
                Result+= sim+(chr(PrefixLength+14)+chr(CodePos-PrefixPos+14))+sim
                CodePos = CodePos + PrefixLength
        return Result
    l_orig=[s_orig[i:i + n] for i in range(0, len(s_orig), n)]
    s_cod = ""
    for s in l_orig:
        Buffer = s
        if Buffer == '':
            break
        s_cod += codeBufferLZ77(Buffer)
    return s_cod
def ILZ77(s_orig, sim):
    DecodedText = ''
    Pos = 0
    i=0
    while i <(len(s_orig)):
        if (s_orig[i] == sim)and(s_orig[i+3] == sim):
            PrefixLength=ord(s_orig[i+1])-14
            Shift =ord(s_orig[i+2])-14
            PrefixPos = Pos-Shift
            DecodedText += DecodedText[PrefixPos:(PrefixPos+PrefixLength)]
            Pos = Pos + PrefixLength
            i+=4
        else:
            DecodedText = DecodedText+s_orig[i]
            Pos = Pos + 1
            i+=1
    return DecodedText
def LZ78(data,sim, max_dict_size: int = 4096,):
    dictionary = {"": 0}
    compressed = []
    prefix = ""
    index = 1
    for char in data:
        new_prefix = prefix + char
        if new_prefix in dictionary:
            prefix = new_prefix
        else:
            compressed.append(chr(dictionary[prefix]+14)+sim+char)
            if len(dictionary) < max_dict_size:
                dictionary[new_prefix] = index
                index += 1
            prefix = ""
    if prefix:
        compressed.append(chr(dictionary[prefix]+14)+sim)
    return "".join(compressed)
def ILZ78(s_cod,sim):
    dictionary = {0: ""}
    decompressed = []
    index = 1
    for entry in range(len(s_cod)):
        if s_cod[entry] == sim:
            prefix_index= ord(s_cod[entry-1])-14
            if entry+1 == len(s_cod):
                char=""
            else:
                char=str(s_cod[entry+1])
            decompressed.append(dictionary[prefix_index] + char)
            dictionary[index] = (dictionary[prefix_index] + char)
            index += 1
    return "".join(decompressed)





def read_byte_str(name, len_sim:int):
    sim_str=""
    with open(name, "rb") as f:
        if len_sim==8:  
            while byte := f.read(1):
                sim_str+=(chr(int(f'{ord(byte):08b}',2)))
        else:
            buf_str=""
            while byte := f.read(1):
                buf_str+="".join(f'{ord(byte):08b}')
                if len(buf_str)>len_sim: 
                    sim_str+=chr(int(buf_str[0:len_sim],2))
                    buf_str=buf_str[len_sim:]
            while len(buf_str)>0:
                if len(buf_str)>len_sim: 
                    sim_str+=chr(int(buf_str[0:len_sim],2))
                    buf_str=buf_str[len_sim:]
                else: 
                    buf_str=buf_str+'0'*(len_sim-len(buf_str))
                    sim_str+=chr(int(buf_str[0:len_sim],2))
                    buf_str=''

        return sim_str
def read_byte_str_bin(name, len_sim:int):
    sim_str=""
    with open(name, "rb") as f:
        if len_sim==8:  
            while byte := f.read(1):
                #print(byte,"-->",f'{ord(byte):08b}')
                sim_str+=f'{ord(byte):08b}'
        else:
            buf_str=""
            while byte := f.read(1):
                buf_str+="".join(f'{ord(byte):08b}')
                if len(buf_str)>len_sim: 
                    sim_str+=buf_str[0:len_sim]
                    buf_str=buf_str[len_sim:]
            while len(buf_str)>0:
                if len(buf_str)>len_sim: 
                    sim_str+=buf_str[0:len_sim]
                    buf_str=buf_str[len_sim:]
                else: 
                    buf_str=buf_str+'0'*(len_sim-len(buf_str))
                    sim_str+=buf_str[0:len_sim]
                    buf_str=''
        return sim_str
def decode_and_read_from_file(filename):
    decoded_string = ""
    with open(filename, 'rb') as file:
        byte = file.read(1)
        while byte:
            byte_value = ord(byte)
            if byte_value < 255:
                decoded_string += chr(byte_value)
                byte = file.read(1)
            else:
                low_byte=0
                while byte_value == 255:
                    low_byte +=  255
                    byte_value=ord(file.read(1))
                decoded_string += chr(byte_value+low_byte)
                byte = file.read(1)
    return decoded_string
def cod_in_byte_str(s):
    print("utf-8cod -completed")
    s_orig=[]  
    for c in s:
        s_orig.append((ord(c)).to_bytes(1, 'big'))
    return s_orig
def print_str_in_file(s,name):
    s=cod_in_byte_str(s)
    with open(name,"wb") as f:
        for l in s:
            f.write(l)
def encode_and_write_to_file(input_string, filename):
    with open(filename, 'wb') as file:
        for char in input_string:
            code = ord(char)
            if code < 255:
                file.write(bytes([code]))
            else:
                while code >=255:
                    file.write(bytes([255]))
                    code-=255
                file.write(bytes([code]))




test=compres_BWT_RLE()
test.compres("test.exe","encoded.txt")
test.decompres("encoded.txt","decompres.exe")
s=decode_and_read_from_file("test.exe")
s_2=decode_and_read_from_file("decompres.exe")
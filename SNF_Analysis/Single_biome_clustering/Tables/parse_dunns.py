file=open("Dunn-test.txt",'r')
file1=open("parsed_results.txt",'w')
l=file.readlines()
sel=[]
for i in l:
    if i[:3]=='[1]':
        sel.append(i)
    elif i[26:-1]=='*':
        sel.append(i[:5])
    else:
        pass
    #file1.write(i)

for i in sel:
    file1.write(str(i)+str("\n"))

file1.close()

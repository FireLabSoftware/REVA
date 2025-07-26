from VSG_Module import *
from sys import argv
from os import getcwd
chrcolor1={'chrI':'red','chrII':'cyan','chrIII':'magenta','chrIV':'blue',
           'chrV':'green','chrX':'black','chrM':'orange','Repetitive':'gray',
           'OP50_Imputed':'yellow','phiX174':'brown'}
try:
    F5n=argv[2]
except:
    from Tkinter import Tk
    root=Tk()
    root.attributes("-topmost", True)
    from tkFileDialog import askopenfilename
    root.withdraw()
    F5n=askopenfilename(title='Incidence Summary File to Open',initialdir=getcwd(), filetypes=[("Text files","IncidenceSummary*.tdv")])
    root.quit()
smoothR1=2
def smootharray(A,r):
    n=len(A)
    A2=[]
    for i in range(n):
        s1=max(0,i-r)
        s2=min(i+r+1,n)
        A2.append(sum(A[s1:s2])/float(s2-s1))
    return A2
    
        


MaxSeparation1=350
vset(bg=white)
Tn5DupMax1=0
F5=open(F5n,mode='rU')
Mnemonic1=F5n.split('_')[2]
NameA1=F5.next().strip().split('\t')[2:]
nN1=len(NameA1)
SeparationArray1=False
for L0 in F5:
    L1=L0.split('\t')
    if L1[0].lower()=='separation':
        L2=map(int,L1[1:])
        if not(SeparationArray1):
            Tn5DupMax1=-L2[0]
            SeparationArray1=[[0]*(Tn5DupMax1+MaxSeparation1+1) for x in range(nN1)]
        if L2[0]<MaxSeparation1+1:
            for i in xrange(nN1):
                SeparationArray1[i][L2[0]+Tn5DupMax1]=L2[i+1]
SeparationArray2=[]
for A1 in SeparationArray1:
    SeparationArray2.append(smootharray(A1,smoothR1))
plotrange1=range(MaxSeparation1+Tn5DupMax1)

for c11 in range(len(NameA1)):
    pr1=max(sum(SeparationArray2[c11]),1)  ## total paired reads
    n1=0
    if pr1<10:
        continue
    if NameA1[c11] in chrcolor1:
        mycolor=chrcolor1[NameA1[c11]]
    else:
        mycolor='e3343'+NameA1[c11]  ## essentially random color for each chromosome 
    for i1 in plotrange1:
        u1=i1-Tn5DupMax1
        ri1=(1.0+SeparationArray2[c11][i1])/pr1           ## normalized tagmentation length incidence
        newxc=(u1+n1)*30; newyc=ri1*100000    
        vcircle(xc=newxc,yc=log(newyc)*1000,r=30,
              fill=mycolor,strokewidth=0,stroke=none,
              xg=u1,yg=ri1,gydom="y2",
              colorkey=NameA1[c11])
        if i1>0:
            pass
            vline(x1=newxc,y1=log(newyc)*1000,x2=oldxc,y2=log(oldyc)*1000, strokewidth=10,stroke=mycolor)
        oldxc=newxc+0; oldyc=newyc+0
vgrid(gtitle="Fragment Length Plot: "+Mnemonic1
      ,gxlabel="Tagmented Length"
      ,gylabel="Relative Incidence"
      ,gylog=True
      ,gxmajor=True
      ,gymajor=True)
vline(x1=(n1)*30,x2=(n1)*30,
      y1=-700,y2=-200,
      strokewidth=50,stroke=cyan,fill=none,
      colorkey='Blunt Junction')
vline(x1=(n1)*30+100,x2=(n1)*30,
      y1=-400,y2=-200,
      strokewidth=50,stroke=cyan)
vline(x1=(n1)*30-100,x2=(n1)*30,
      y1=-400,y2=-200,
      strokewidth=50,stroke=cyan)
vline(x1=(n1+8)*30,x2=(n1+8)*30,
      y1=-700,y2=-200,
      strokewidth=50,stroke=orange,fill=none,
      colorkey='Full Length + 9b overhang')
vline(x1=(n1+8)*30+100,x2=(n1+8)*30,
      y1=-400,y2=-200,
      strokewidth=50,stroke=orange)
vline(x1=(n1+8)*30-100,x2=(n1+8)*30,
      y1=-400,y2=-200,
      strokewidth=50,stroke=orange)
vcolorkey()
vdisplay(Mnemonic1.split('.')[0]+'A_FragHistogram.svg')

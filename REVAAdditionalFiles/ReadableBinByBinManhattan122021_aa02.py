## A script that uses VSG_Module to prpduce a Manhattan plot of "Experimental" versus "Total_Genomic" (e.g. eccDNA/Total_Genomic) with each normalized expectations from genome k-mer counts on 
MyVariables = {'ufn':'UnEnrichedFileName',
                'efn':'EnrichedFileName',
                'bs':'BinSize',
                'mukpb':'MinUniqueKmersPerBin',
                'muspb':'MinUniqueStartsPerBin', 
                'rkc':'ReferenceKmerColumn', 
                'dskcc':'DataSenseKmerCountColumn', 
                'dakcc':'DataAntisenseKmerCountColumn', 
                'dskcc':'DataSenseKmerCountColumn', 
                'dakcc':'DataAntisenseKmerCountColumn', 
                'dsssc':'DataSenseStartSpeciesColumn',  
                'dsssc':'DataAntisenseStartSpeciesColumn',  
                'hbtpr':'HorizontalBinToPixelRatio', 
                'vcsp':'VerticalChromsomeSeparationPixels',
                'fetpr':'FoldEnrichmentToPixelRatio',
                'kcrv':'KmerCountRegularizationValue',
                'srv':'SpeciesRegularizationValue'} 
for mv in list(MyVariables.keys()):
    MyVariables[MyVariables[mv].lower()]=mv

UnEnrichedFileName = 'BinByBinReadCountSummary_REVA_CAT3_S3_L001_R2_001_AllHuman38_NO_Alts_D_11_11_21_T_08_16_50.tdv'
EnrichedFileName = 'BinByBinReadCountSummary_REVA_CAT3exoII_S19_L001_R2_001_AllHuman38_NO_Alts_D_11_11_21_T_07_56_55.tdv'
BinSize = 100000
MinUniqueKmersPerBin = 1000
MinUniqueStartsPerBin = 10
ReferenceKmerColumn = 3
DataSenseKmerCountColumn = 8
DataAntisenseKmerCountColumn = 10
DataSenseStartSpeciesColumn = 11
DataAntisenseStartSpeciesColumn = 13
HorizontalBinToPixelRatio = 2
VerticalChromsomeSeparationPixels = 250
FoldEnrichmentToPixelRatio = 5
KmerCountRegularizationValue = 100
SpeciesRegularizationValue = 5
from sys import argv
GraphLegends = ['Program = ' + argv[0]]
for a1 in argv[1:]:
    VarMnemonic = a1.split('=')[0].strip()
    VarName = MyVariables[VarMnemonic.lower()]
    VarValue = a1.split('=')[1].strip()
    if VarValue.isnumeric():
        locals()[VarName] = int(VarValue)
    else:
        locals()[VarName] = VarValue
    GraphLegends.append(VarName+' = '+VarValue)
        

def tdvToDict(Fn1):
## This routine will take a tab-delimited-value file and return a dictionary where the keys are gene names parsed from each line and the values are lists of columns
    F = open(Fn1,mode='rt').readlines()
    D = {} ## This will be the output dictionary
    for L0 in F:
        if L0.startswith('chr'):
            L1 = L0.strip().split('\t') ## L1 will be a list of entries derived from the tab-delimited line
            g = L1[0]+':'+L1[1] ## This line parses the first two items in a bin definition (chromsome and start) into a unique identifier
            L2 = [] ## This will be a list where integer values are represented by integer objects
            for x in L1:
                if x.isnumeric():
                    L2.append(int(x)) 
                else:
                    L2.append(x) 
            D[g] = L2
    return D

from VSG_ModuleDD import *
vset(bg=white)

D1 = tdvToDict(UnEnrichedFileName) ## Data Dictionary for Unenriched data.  Keys are BinID strings, values are lists of columns
D2 = tdvToDict(EnrichedFileName) ## Data Dictionary for Enriched data
tk1 = sum([D1[x][ReferenceKmerColumn] for x in D1]) ## Total reference genome unique Kmer count
tcG1 = sum([D1[x][DataSenseKmerCountColumn] for x in D1])+sum([D1[x][DataAntisenseKmerCountColumn] for x in D1]) ## Total Kmer matches in Unenriched data
tcE1 = sum([D2[x][DataSenseKmerCountColumn] for x in D2])+sum([D2[x][DataAntisenseKmerCountColumn] for x in D2]) ## Total Kmer matches in Enriched data
tsG1 = sum([D1[x][DataSenseStartSpeciesColumn] for x in D1])+sum([D1[x][DataAntisenseStartSpeciesColumn] for x in D1]) ## How many different unique start sites are seen in Unenriched data
tsE1 = sum([D2[x][DataSenseStartSpeciesColumn] for x in D2])+sum([D2[x][DataAntisenseStartSpeciesColumn] for x in D2]) ## How many different unique start sites are seen in Enriched data
reg_tk1 = KmerCountRegularizationValue + tk1
reg_tcG1 = KmerCountRegularizationValue + tcG1
reg_tcE1 = KmerCountRegularizationValue + tcE1
reg_tsG1 = SpeciesRegularizationValue + tsG1
reg_tsE1 = SpeciesRegularizationValue + tsE1

## Calculate Medians Coverage
eEList1 = []
for g1 in D1:
    sG1 = D1[g1][DataSenseStartSpeciesColumn]+D1[g1][DataAntisenseStartSpeciesColumn]
    sE1 = D2[g1][DataSenseStartSpeciesColumn]+D2[g1][DataAntisenseStartSpeciesColumn]
    if sG1>=MinUniqueStartsPerBin:
        eEList1.append(sE1/sG1)
Med_eE1 = sorted(eEList1)[len(eEList1)//2]  ## Median ratio of reads 

SizeByC1 = {}
for g1 in D1:
    k1 = D1[g1][ReferenceKmerColumn]
    cG1 = D1[g1][DataSenseKmerCountColumn]+D1[g1][DataAntisenseKmerCountColumn]
    cE1 = D2[g1][DataSenseKmerCountColumn]+D2[g1][DataAntisenseKmerCountColumn]
    sG1 = D1[g1][DataSenseStartSpeciesColumn]+D1[g1][DataAntisenseStartSpeciesColumn]
    sE1 = D2[g1][DataSenseStartSpeciesColumn]+D2[g1][DataAntisenseStartSpeciesColumn]
    reg_k1 = KmerCountRegularizationValue + k1
    reg_cG1 = KmerCountRegularizationValue + cG1
    reg_cE1 = KmerCountRegularizationValue + cE1
    reg_sG1 = SpeciesRegularizationValue + sG1
    reg_sE1 = SpeciesRegularizationValue + sE1
    if k1<MinUniqueKmersPerBin or sG1<MinUniqueStartsPerBin: ## This is an important statement that filters out features where there are too few unique k-mers (k1<100) or the number of different start positions in genomic DNA is small (<10) 
        continue
    Amp1 = (reg_cG1/reg_k1)/(reg_tcG1/reg_tk1) ## This is an estimate of KMers_Counts_Observed/KMers_Counts_Expected for the Unenriched sample
    ##eE1 = (reg_cE1/reg_tcE1)/(reg_cG1/reg_tcG1) ## k-mer-count based enrichment ratio in enriched sample: ratio of obs/expected in enriched, divided by same ratio on gDNA
    eE1 = (sE1/sG1)/Med_eE1 ##(reg_sG1/reg_tsG1) ## unique-read-species-count based enrichment in enrichment ratio in enriched sample, ratio of obs/expected in ecc, divided by same ratio on gDNA
    chr1 = D1[g1][0]
    c1 =  chr1[3:] ## the chromsome number (zero for M, 23 for X, 24 for Y)
    if c1.startswith('EBV') or c1.startswith('HBV') or c1.startswith('Y') or c1.startswith('M'): continue
    if c1=='X':
        chrnum = 23
    elif c1=='Y':
        chrnum = 24
    else:
        chrnum = int(c1)
    y0 = -chrnum*VerticalChromsomeSeparationPixels 
    SizeByC1[chr1] = (y0,D1[g1][1]+D1[g1][2])
    p0 = int(g1.split(':')[1]) ## lateral position in chromosome of first base
    p00 = max(1,p0-2*BinSize)
    p1 = D1[g1][1]+D1[g1][2]-1
    p11 = p1+2*BinSize
    color1 = Amp1/(Amp1+1.0)
    xc1 = p0*HorizontalBinToPixelRatio/BinSize
    yc1 = min(VerticalChromsomeSeparationPixels*0.85,eE1*FoldEnrichmentToPixelRatio)
    xlink1 = "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position="+chr1+"%3A"+str(p00)+"%2D"+str(p11)+"&highlight="+chr1+"%3A"+str(p0)+"%2D"+str(p1)
    vcircle(xc=xc1,
            yc=y0+yc1,
            fill=color1,
            stroke=black,
            strokewidth=2,
            r=2*HorizontalBinToPixelRatio,
            xlink=xlink1,
            strokeopacity=0.3)
    
for chr1 in SizeByC1:
    if chr1=='chrM': continue
    vtext(text=chr1,x1=SizeByC1[chr1][1]*HorizontalBinToPixelRatio/BinSize+5*HorizontalBinToPixelRatio, yc=SizeByC1[chr1][0]+VerticalChromsomeSeparationPixels/2, font="DejaVuSerif 54 Bold", fill=black)
    vrect(x1=0, x2=SizeByC1[chr1][1]*HorizontalBinToPixelRatio/BinSize, y1=SizeByC1[chr1][0], y2=SizeByC1[chr1][0]+VerticalChromsomeSeparationPixels, stroke=black, fill=none, strokewidth=2, priority=-10)
HyperG1 = lambda x:vcolor(x/(1+x))
vlegend(font="DejaVuSerif 36 Bold")
vlegend(text='\r'.join(GraphLegends), font="DejaVuSerif 36 Bold")
vcolorkey(spectrumtitle='Amplification In Unenriched Sample', colorvalue=HyperG1, mincolorindex=0.099, maxcolorindex=10, logmode=True)
vdisplay('BinByBinManhattanPlot001.svg')

__author__ = 'lucas'
import re
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
motivos = open('teste.mtv','r').read()
m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
for linha in motivos.split('\n'):
    m.add_instance(Seq(linha,m.alphabet))
    print linha

print m.background
#m.weblogo('logo.png')




motivos_inicio_0_7=[]
arquivo = open('/home/lucas/PycharmProjects/hmm_cirRNA/original_sem_deletar_os_com_gap_introns_ciRNAs.txt','r').read()
for seq in arquivo.split('>'):
    if len(seq)>0:
        header= seq.split('\n')[0]
        position= header.split('_')[6].split('-')
        strand= header.split('_')[7]
        branch= int(re.sub(';','',header.split('_')[8]))
        #print branch,header
        #print strand,position
    for seq_bases in re.split('hg19.*\n' ,seq):
        sequencia_completa = re.sub('\r|\n','',seq_bases)

    #print sequencia_completa

    if sequencia_completa:

        for position,sequence in  m.search_pwm(Seq(sequencia_completa,m.alphabet)):
            print position
            print sequencia_completa[position:position+11]
            #print 'teste'
           # print (position,sequence.tostring())
            #if position>0:
            #print sequencia_completa[position:position+11]
            #else:
            #    print sequencia_completa[position:position+11]
        #print sequencia_completa[0:7]
        #motivos_inicio_0_7.append(sequencia_completa[0:7])

        #print int(position[1]),branch
        #branch_bp=int(position[1])-branch
        #print branch_bp
        #bp_exact=len(sequencia_completa)-branch_bp
        #print branch_bp
        #print sequencia_completa[bp_exact-1]
        #print strand,position,branch
        #print header
        #print sequencia_completa
        #if strand=='+':
        #    branch_bp=int(position[1])-int(branch)
            #print branch_bp
        #    bp_exact=len(sequencia_completa)-branch_bp
        #    print branch_bp
         #   print '>seq'+str(branch)
         #   print sequencia_completa[bp_exact-13:bp_exact]
        #else:
        #    branch_bp=branch-int(position[0])
            #print branch_bp
        #    bp_exact=len(sequencia_completa)-branch_bp
            #print bp_exact
        #    print '>seq'+str(branch)
        #    print sequencia_completa[bp_exact-13:bp_exact]
            #print branch_bp
        #    print branch,sequencia_completa[bp_exact-1]
        #    branch_bp= int(position[0])-branch
        #    print sequencia_completa[branch_bp]
        #else:
        #    branch_bp= int(position[1])-branch
        #    print sequencia_completa[branch_bp]

        #for bp in re.findall('\w',sequencia_completa):
        #    print bp
        #


        #print sequencia_completa
#print motivos_inicio_0_7
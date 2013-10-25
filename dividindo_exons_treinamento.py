import re
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
arquivo = open('introns_ciRNAs.txt','r').read()
frequencias_5=[]
frequencias_3=[]

separando= re.split('>',arquivo)
for seq in separando:
    if len(seq)>0:
        #print seq
        header=seq.split(r';')[0]
        seq=seq.split(r';')[1]
        seq=re.sub('\n|\r','',seq)

        #print seq
        print '>'+header
        calculo = header.split('_')[6].split('-')

        calculo.append(header.split('_')[8])
        #print calculo
        if header.split('_')[7]=='+':
            inicio = (int(calculo[1]))-int(calculo[2])
        else:
            inicio = (int(calculo[0]))-int(calculo[2])
        inicio= abs(inicio)
        inicio=len(seq)-inicio
        print seq[inicio-11:inicio]
        m.add_instance(Seq(seq[inicio-11:inicio],m.alphabet))
for x in  m.pwm():
    print x
m.weblogo('teste.png')
    #print seq

#    for corte in seq.split('\n'):
#        #print corte
#        corte=re.sub('\r|\n','',corte)
#        if len(corte)>1:
#            #print '<'+corte+'>'
#            #pass
#            if  not re.match(r'hg',corte):
#                adicionar=adicionar+corte
#
#
#
#                print 'entrou', corte
#                #print seq[0:6],'3 linha'
#                #print seq[len(seq)-11:len(seq)],'5 linha'
#
#                #m.add_instance(Seq(seq[0:7],m.alphabet))
#
#            #else:
#                #print corte
#                #pass
#        print adicionar
#        print adicionar[0:6],'3 linha'
#        print adicionar[len(seq)-11:len(seq)],'5 linha'
#        adicionar=adicionar[len(adicionar)-11:adicionar(adicionar)]
#
#
#
#for x in  m.pwm():
#   # print x
#
#    pass

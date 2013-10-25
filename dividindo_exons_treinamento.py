import re
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
arquivo = open('introns_ciRNAs.txt','r').read()
frequencias_5=[]
frequencias_3=[]
treinamento=[]
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
        entrada_vetor=[]
        ponta=0
        if header.split('_')[7]=='+':
            inicio = (int(calculo[1]))-int(calculo[2])
            #ponta=int(calculo[1])
        else:
            inicio = (int(calculo[0]))-int(calculo[2])
            #ponta=int(calculo[0])
        inicio= abs(inicio)
        inicio=len(seq)-inicio
        #print seq[inicio-11:inicio]
        #print seq[0:6]
        m.add_instance(Seq(seq[inicio-11:inicio],m.alphabet))
        #m.add_instance(Seq(seq[0:6],m.alphabet))
        onze=11
        for i in range(len(seq)):

            '''
            capturando a ponta ate  o branch
            '''
            if i > inicio:
                print ponta,'ponta',seq[i],i
                entrada_vetor.append('antes_branch')
                entrada_vetor.append(seq[i])
            '''
            capturando branch
            '''
            if i == inicio:
                print ponta,'BRANCH',seq[i],i
                entrada_vetor.append('BRANCH')
                entrada_vetor.append(seq[i])

            if i < inicio and i > (inicio-11):
                print ponta,'branch_'+str(onze),seq[i],i
                entrada_vetor.append('branch_'+str(onze))
                entrada_vetor.append(seq[i])
                onze=onze-1
            if i <(inicio-11) and i>7:
                print ponta,'corpo',seq[i],i
                entrada_vetor.append('corpo')
                entrada_vetor.append(seq[i])
            if i <= 7:
                print ponta,'FINAL',seq[i],i
                entrada_vetor.append('Final')
                entrada_vetor.append(seq[i])
        print entrada_vetor
for x in  m.pwm():
    print x
m.weblogo('teste2.png')
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

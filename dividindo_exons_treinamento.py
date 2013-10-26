import re
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import simplehmm
m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
arquivo = open('introns_ciRNAs.txt','r').read()
frequencias_5=[]
frequencias_3=[]
treinamento=[]
entrada_treinamento=[]
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
        final=1
        '''CHECAR OS ESTADOS FINAIS FINAL8'''
        estados=['antes_branch','BRANCH','branch_11','branch_10','branch_9','branch_8','branch_7','branch_6','branch_5','branch_4','branch_3','branch_2','corpo','Final1','Final2','Final3','Final4','Final5','Final6','Final7','Final8']
        emissoes=['A','C','T','G']
        treinamento=[]
        for i in range(len(seq)):
            debug=seq[i]
            entrada_vetor=[]
            '''
            capturando a ponta ate  o branch
            '''
            if i > inicio:
                #print ponta,'ponta',seq[i],i

                entrada_vetor.append(seq[i])
                entrada_vetor.append('antes_branch')
                treinamento.append(entrada_vetor)
            '''
            capturando branch
            '''
            if i == inicio:
                #print ponta,'BRANCH',seq[i],i

                entrada_vetor.append(seq[i])
                entrada_vetor.append('BRANCH')
                treinamento.append(entrada_vetor)
            '''
            capturando_11 apos branch
            '''
            if i < inicio and i > (inicio-11):
                #print ponta,'branch_'+str(onze),seq[i],i

                entrada_vetor.append(seq[i])
                entrada_vetor.append('branch_'+str(onze))
                onze=onze-1
                treinamento.append(entrada_vetor)
            if i <(inicio-11) and i>7:
                #print ponta,'corpo',seq[i],i

                entrada_vetor.append(seq[i])
                entrada_vetor.append('corpo')
                treinamento.append(entrada_vetor)
            if i <= 7:
               # print ponta,'FINAL'+str(final),seq[i],i

                entrada_vetor.append(seq[i])
                entrada_vetor.append('Final'+str(final))
                final+=1
                treinamento.append(entrada_vetor)
            print entrada_vetor
    if len(treinamento):
        entrada_treinamento.append(treinamento)

for y in entrada_treinamento:
    print y





#for y in entrada_treinamento:
#    print y
cirRNA_hmm= simplehmm.hmm('cirRNA_hmm_primeiro_teste',estados,emissoes)
#teste=simplehmm.hmm()
cirRNA_hmm.train(entrada_treinamento, smoothing='absdiscount')
print cirRNA_hmm.print_hmm()
#for imprimir in treinamento:
#    print imprimir
#for x in  m.pwm():
#    print x
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

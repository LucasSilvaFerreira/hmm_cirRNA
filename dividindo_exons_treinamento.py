import re
from Bio import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import simplehmm
from random import random, shuffle

m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
arquivo = open('introns_ciRNAs.txt','r').read().upper()
#arquivo = open('query.txt','r').read().upper()
frequencias_5=[]
frequencias_3=[]
treinamento=[]
entrada_treinamento=[]
query=[]
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
        inicio-= 1 # the minus one is put because the vector starts in zero
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
        query.append(list(seq))
        contador=1
        for i in range(len(seq)):

            debug=seq[i]
            entrada_vetor=[]
            '''
            capturando a ponta ate  o branch



            '''





            if i > inicio:
                #print ponta,'ponta',seq[i],i


                entrada_vetor.append('antes_branch')
                entrada_vetor.append(seq[i])

                treinamento.append(entrada_vetor)
            '''
            capturando branch
            '''
            if i == inicio:
                print ponta,'BRANCH',seq[i],i

                print contador
                entrada_vetor.append('BRANCH')
                entrada_vetor.append(seq[i])
                treinamento.append(entrada_vetor)
            '''
            capturando_11 apos branch
            '''
            if i < inicio and i > (inicio-11):
                #print ponta,'branch_'+str(onze),seq[i],i


                entrada_vetor.append('branch_'+str(onze))
                entrada_vetor.append(seq[i])
                onze=onze-1
                treinamento.append(entrada_vetor)
            if i <(inicio-11) and i>7:
                #print ponta,'corpo',seq[i],i


                entrada_vetor.append('corpo')
                entrada_vetor.append(seq[i])

                treinamento.append(entrada_vetor)
            if i <= 7:
               # print ponta,'FINAL'+str(final),seq[i],i


                entrada_vetor.append('Final'+str(final))
                entrada_vetor.append(seq[i])
                final+=1
                treinamento.append(entrada_vetor)
            #print entrada_vetor
            contador+=1
    if len(treinamento):
        entrada_treinamento.append(treinamento)

'''TREINAMENTO

'''


#cirRNA_hmm= simplehmm.hmm('cirRNA_hmm_primeiro_teste',estados,emissoes)
#cirRNA_hmm.train(entrada_treinamento, smoothing='absdiscount')
#cirRNA_hmm.save_hmm("circ.hmm")
#print cirRNA_hmm.print_hmm()



'''VALIDACAO'''

cirRNA_hmm=simplehmm.hmm('circ_rna primeiros testes',['dummy'], ['dummy'])
cirRNA_hmm.load_hmm('circ.hmm')
#print query[1]
#print cirRNA_hmm.print_hmm()


print len(query[0]),'size'
for query_out in query:
    print query_out



print cirRNA_hmm.viterbi(query[0])[1]


print "------------------------------------------"
print cirRNA_hmm.viterbi(query[0])[0]
#print cirRNA_hmm.viterbi(query[0])[0].index('BRANCH'),'<-BRANCH->',query[0][cirRNA_hmm.viterbi(query[0])[0].index('BRANCH')]
for i in range(1,10):
    shuffle(query[0])
    print cirRNA_hmm.viterbi(query[0])[1],"---"

#print len(query),'query'


#for testar in query:
#    for embaralhar in range(1,2):
#        shuffle(testar)
#        try:
#            #print cirRNA_hmm.viterbi(testar)[0]
#            #print cirRNA_hmm.viterbi(testar)[0].index('BRANCH'),'<-BRANCH->',query[0][cirRNA_hmm.viterbi(query[0])[0].index('BRANCH')]
#            print cirRNA_hmm.viterbi(testar)[1]
#            #print cirRNA_hmm.viterbi(testar)[0]
#        except:
#            print ('probabilidade baixa')
#

#print query[0]

print query[0]

#for imprimir in treinamento:
#    print imprimir
#for x in  m.pwm():
#    print x
#m.weblogo('teste2.png')



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

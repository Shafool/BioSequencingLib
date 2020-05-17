# -*- coding: utf-8 -*-
"""
Created on Tue May  5 20:44:22 2020

@author: Acer
"""

#Clase de secuencia biologica
class Myseq:
    
   
        
    def __init__(self, seq, tipoSecuencia="DNA"):
        self.seq=seq
        self.tipoSecuencia=tipoSecuencia
    def imprimir_secuencia(self):
        print("Secuencia:",self.seq)
    def getbiotiposeq(self):
        return self.tipoSecuencia
    def MostrarInfoSeq(self):
        # el tipo de secuencia se actualiza a un biotipo inválido 
        # por alteración directa del atributo
        s1.tipoSecuencia  = "time series"
        
        #Alternativa segura: metodo para validar el cambio
    def setbiotiposeq(self, bt):
        biotipo = bt.upper()
        if biotipo == "DNA" or biotipo == "RNA" or biotipo == "PROTEINA":
             print("Biotipo: ",biotipo)
             self.tipoSecuencia=biotipo
            
               
        else:
            print("Secuencia no biologica!")
        
#s1.setbiotiposeq("time series")
#s1.setbiotiposeq("dna")
#print(len(s1))
#print(s1[4])
#print(s1[2:5])

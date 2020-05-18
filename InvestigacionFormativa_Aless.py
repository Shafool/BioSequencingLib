# -*- coding: utf-8 -*-
"""
Created on Sat May 16 20:06:48 2020

@author: Acer
"""

import sys,time
#Clase de secuencia biologica
class MySeq:
    
     
    def __len__(self):
        return len(self.seq)  
    
    def __str__(self):
        return self.seq
    
    def __getitem__(self,n):
        return  self.seq[n]
    
    def __init__(self, seq, tipoSecuencia=""):
        self.seq=seq
        self.tipoSecuencia=tipoSecuencia
    def imprimir_secuencia(self):
        print("Secuencia:",self.seq)
    def getbiotiposeq(self):
        return self.tipoSecuencia
    def MostrarInfoSecuencia(self):
        texto = str(s1)
        cadena = list(texto)
        tamanio = len(cadena)
        ListaADN=['A','G','C','T']
        ListaSecuencia=cadena
        tamanio2= len(ListaADN)
        secuencia = self.seq
        
        if ('R' in secuencia) or ('N' in secuencia) or ('D' in secuencia) or ('E' in secuencia) or ('Q' in secuencia) or ('H' in secuencia) or ('I' in secuencia) or ('L' in secuencia) or ('K' in secuencia) or ('M' in secuencia) or ('F' in secuencia) or ('P' in secuencia) or ('S' in secuencia) or ('W' in secuencia) or ('Y' in secuencia) or ('V' in secuencia) :
            s1.setbiotiposeq("PROTEINA")
            
                                                                   
        
        if 'U' in secuencia:
            s1.setbiotiposeq("RNA")
         
        
        
        if 'A' and 'G'and 'C'and'T' in ListaSecuencia and 'U' not in ListaSecuencia:
            for i in range (0,tamanio-1):
                for j in range(0, tamanio2-1):
                    if ListaSecuencia[i]==ListaADN[j]:
                        i=i+1
                    else:
                        
                        if j==tamanio2-1:
                                while j!=0:
                                    
                                   j=j-1
                         
                                
                                
                        else:
                            if (( (j-1)!=-1) and j+1==4):
                                j=j-1
                            else:
                                j=j+1
                            
            s1.setbiotiposeq("DNA")        
                        
       
            
        
        
    def ContarOcurrencias(self, buscar_seq):
        return self.seq.count(buscar_seq)
       
    def setbiotiposeq(self, bt):
         biotipo = bt.upper()
         if biotipo == "DNA" or biotipo == "RNA" or biotipo == "PROTEINA":
             print("Biotipo: ",biotipo)
             self.tipoSecuencia=biotipo
             
             texto = str(s1)
             cadena = list(texto)
             tamanio = len(cadena)
             
             #print ("Cadena de ADN valida: ",cadena)          
              #"""RNA"""
             if biotipo == "RNA":
                 
                  print("Validando componentes...")
                  print ("Cadena: ",cadena)
                  nroA=cadena.count("A")
                  nroC=cadena.count("C")
                  nroG=cadena.count("G")
                  nroT=cadena.count("T")
                  nroU=cadena.count("U")
                  
                  print(nroA,"valores de A Identificadas en la secuencia")
                         
                
                         
                  print(nroC,"valores de C Identificadas en la secuencia")
                         
                     
                         
                  print(nroG,"valores de G Identificadas en la secuencia")
                          
                     
                         
                  print(nroT,"valores de T Identificadas en la secuencia")
                  
                  if 'U' in cadena:
                      
                      print(nroU,"valores de U Identificadas en la secuencia")
                  
                  if 'U' in cadena and 'T' not in cadena:
                      print("Secuencia valida: ", cadena)
                      print("Finalizando programa")
                      sys.exit()
                      
                  
                  
                  while True:
                      
                      try:
                          if ('U' in cadena and 'T' in cadena) or ('U' not in cadena) or ('T' not in cadena):
                              
                              rpta = input("No se encontro U en la secuencia o faltan valores por convertir, desea convertir T en U? (Si / No):")
                         
                      except ValueError:
                          print("Ingrese una respuesta valida")
                          continue
                      else:
                          break
                          
                  
                      
                  if (rpta == "Si"):
                      print("Convirtiendo...")
                      for n in range(0,tamanio):
                          if cadena[n]=="T":
                              cadena[n]="U"
                          time.sleep(1)

                          print (cadena[n])
                          time.sleep(1)
                          print ("Nueva cadena: ",cadena)
                      else:
                          print("Conversion fizalizada")
                          s1.seq=''.join(cadena)
                          s1.setbiotiposeq("RNA")
                          
                          
                          
                          print("Secuencia de RNA:",cadena)
                          
                          
                          
                  if (rpta=="No"):
                      print("No se encuentra U o faltan valores por convertir, No es RNA")
                      sys.exit()

                          
                          
                  
                 
                          
                  
                  
              #"""PROTEINA"""
             if biotipo == "PROTEINA":
                 print("Futuras operaciones con proteinas...")
                 sys.exit()
                 
                  
                      
                  
                  
                  
             #"""DNA"""
             if biotipo == "DNA":
                 
                 print("Validando componentes...")
                 print ("Cadena: ",cadena)
                 nroA=cadena.count("A")
                 nroC=cadena.count("C")
                 nroG=cadena.count("G")
                 nroT=cadena.count("T")
                 for n in range(0,tamanio):
                     
                     if 'A' in cadena and 'C' in cadena and 'G' in cadena and 'T' in cadena :
                         
                         print("Cadena valida!")
                         break     
                     else:
                         print("Cadena no valida!")
                         break 
                
                   
                 print(nroA,"valores de A Identificadas en la secuencia")
                         
                
                         
                 print(nroC,"valores de C Identificadas en la secuencia")
                         
                     
                         
                 print(nroG,"valores de G Identificadas en la secuencia")
                         
                     
                         
                 print(nroT,"valores de T Identificadas en la secuencia")
                    
                     
            
             
                 
                 
         else:
            print("Secuencia no biologica!, intente reescribir la secuencia o el tipo!")
        
    
    
 
#Especificando la secuencia
Secuencia = input("Ingrese una secuencia: ")            
s1 = MySeq(Secuencia)  
"""CAATAGCUACTG"""
#LLamando metodos
s1.imprimir_secuencia()
s1.MostrarInfoSecuencia()
#Especificando el biotipo
    
    
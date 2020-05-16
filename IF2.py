# -*- coding: utf-8 -*-
"""
Created on Thu May 14 23:22:53 2020

@author: Acer
"""

from sequences import*


def validarADN(secuencia_ADN):
    
#Valida si la secuencia de ADN es valida. Devuelve True si es valida False si no
    seqm =secuencia_ADN.upper()
    valido = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")
    if valido == len(seqm): return True
    else: return False


seq = input("Ingrese secuencia de ADN: ")
if validarADN(seq):
    print("Secuencia Valida:")
    print("Transcripcion:")
    print("Complemento Inverso:")
    print("Contenido global de GC:")
    print("Translacion directa:")
    print("Todas las proteinas en el ORFS (tamaño decrementado):")
    
else: print("La secuencia de ADN no es valida")
    
    